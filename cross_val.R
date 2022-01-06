library(here)
library(reticulate)

source(paste0(here(), "/load_data.R"))

source_python("learn_metrics.py")

# melting data from Table S7
mdata <- load_mdata()

# convert melting data to parametric melting data
mparams <- param_mdata(mdata) %>% drop_na()

# pairwise interactions from table S2
ints <- load_interactions() %>% mutate(Num_pubs = Num_pubs %>% as.integer())

# Remove each interaction where at least one protein does not have parametric melting data
ints <- ints %>%
  semi_join(mparams, by = c("ProteinA" = "Protein")) %>%
  semi_join(mparams, by = c("ProteinB" = "Protein"))

# Ensure there are no duplicate interactions
ints <- ints %>%
  mutate(Proteins = map2(ProteinA, ProteinB, ~ sort(c(.x, .y)))) %>%
  mutate(ProteinA = map_chr(Proteins, ~ .x[1])) %>%
  mutate(ProteinB = map_chr(Proteins, ~ .x[2])) %>%
  distinct(ProteinA, ProteinB, .keep_all = T) %>%
  select(ProteinA, ProteinB, Num_pubs)

# CORUM complex data from table S8
cdata <- load_cdata()

# Remove each protein with no parametric melting data
cdata <- cdata %>% semi_join(mparams, by = "Protein")

# pairs of proteins that belong to the same complex
cints <- cdata %>%
  nest(Proteins = -Complex) %>%
  mutate(Proteins = map(Proteins, ~ .x$Protein)) %>%
  mutate(Interactions = map(Proteins, ~{
    m <- t( combn(.x, 2) );
    colnames(m) <- c("ProteinA", "ProteinB");
    m %>% as_tibble()
  })) %>%
  select(Interactions, Complex) %>%
  unnest(Interactions)

# Ensure there are no duplicate interactions
cints <- cints %>%
  mutate(Proteins = map2(ProteinA, ProteinB, ~ sort(c(.x, .y)))) %>%
  mutate(ProteinA = map_chr(Proteins, ~ .x[1])) %>%
  mutate(ProteinB = map_chr(Proteins, ~ .x[2])) %>%
  distinct(ProteinA, ProteinB, .keep_all = T) %>%
  select(ProteinA, ProteinB, Complex)

# all_ints: union of ints and cints
all_ints <- bind_rows(ints %>% select(-Num_pubs), cints %>% select(-Complex)) %>%
  distinct(ProteinA, ProteinB)

# Each protein's interactors
interactors <- bind_rows(
  all_ints,
  all_ints %>% select(ProteinA = ProteinB, ProteinB = ProteinA)
) %>%
  nest(Interactors = -ProteinA) %>%
  mutate(Interactors = map(Interactors, ~ .x$ProteinB %>% unique)) %>%
  select(Protein = ProteinA, Interactors)

# cross-validation interactions
cv_ints <- ints %>% filter(Num_pubs >= 2) %>% select(ProteinA, ProteinB)

# Randomly label the interactions so that they can be partitioned into five sets of approximately the same
# size. 
set.seed(27)
cv_ints <- cv_ints %>%
  mutate(Label = rep(1:5, ceiling(nrow(cv_ints) /5)) %>% sample(., size = nrow(cv_ints)))

# convert interactions to training data
make_training_data <- function(ints) {
  
  # pairs of non-interacting proteins
  pairs_ni <- tibble(ProteinA = sample(mparams$Protein, size = nrow(ints), replace = T)) %>%
    left_join(interactors, by = c("ProteinA" = "Protein")) %>%
    mutate(ProteinB = map2_chr(ProteinA, Interactors, ~{
      repeat {
        protB <- sample(mparams$Protein, 1)
        if ( !(protB %in% c(.x,.y)) ) {
          break
        }
      }
      protB;
    })) %>%
    select(-Interactors)
  
  pairs <- bind_rows(ints, pairs_ni)
  
  prots <- c(pairs$ProteinA, pairs$ProteinB) %>% unique
  
  # melting data
  md <- tibble(Protein = prots, Index = 0:(length(prots) - 1), ) %>%
    inner_join(mdata, by = "Protein")
  
  # X
  X <- md %>% select(T40:T64) %>% as.matrix()
  colnames(X) <- NULL
  # Convert to list of vectors
  X <- array_tree(X, margin = 1)
  
  # pairs_indices
  pairs_indices <- pairs %>% 
    inner_join(md %>% select(Protein, Index), by = c("ProteinA" = "Protein")) %>%
    inner_join(md %>% select(Protein, Index), by = c("ProteinB" = "Protein")) %>%
    select(Index.x, Index.y) %>%
    as.matrix()
  colnames(pairs_indices) <- NULL
  pairs_indices <- array_tree(pairs_indices, margin = 1)
  
  # y_pairs
  y_pairs <- c(rep(1L, nrow(ints)), rep(-1L, nrow(pairs_ni)))
  
  return(list(X = X, pairs_indices = pairs_indices, y_pairs = y_pairs))
}

# Use cross-validation to learn a metric for each fold
learn_metrics_cv <- function(cv_ints) {
  
  lbls <- cv_ints$Label %>% unique %>% sort
  itml_metrics <- vector(mode = "list", length = length(lbls))
  itml_neg_metrics <- vector(mode = "list", length = length(lbls))
  mmc_neg_metrics <- vector(mode = "list", length = length(lbls))
  mmc_diag_neg_metrics <- vector(mode = "list", length = length(lbls))
  for (l in lbls) {
    
    tr <- cv_ints %>%
      filter(Label != l) %>%
      select(-Label) %>%
      make_training_data()
    
    print(paste0("Learning ITML metric for label ", l, "..."))
    itml_metrics[[l]] <- learn_itml(tr, min_num = length(tr$y_pairs), use_neg_pairs = F, verbose = T)
    
    print(paste0("Learning ITML-neg metric for label ", l, "..."))
    itml_neg_metrics[[l]] <- learn_itml(tr, min_num = length(tr$y_pairs), use_neg_pairs = T, verbose = T)
    
    print(paste0("Learning MMC-neg metric for label ", l, "..."))
    mmc_neg_metrics[[l]] <- learn_mmc(tr, min_num = length(tr$y_pairs), use_neg_pairs = T, diag = F, initialization = "identity", verbose = T, max_iter = 200L)
    
    print(paste0("Learning MMC-diag-neg metric for label ", l, "..."))
    mmc_diag_neg_metrics[[l]] <- learn_mmc(tr, min_num = length(tr$y_pairs), use_neg_pairs = T, diag = T, initialization = "identity", verbose = T)
  }
  
  return(list(ITML = itml_metrics, ITML_neg = itml_neg_metrics, MMC_neg = mmc_neg_metrics, MMC_diag_neg = mmc_diag_neg_metrics))
}

cv_metrics <- learn_metrics_cv(cv_ints)
write_rds(cv_ints, paste0(here(), "/cv_ints.rds"))
write_rds(cv_metrics, paste0(here(), "/cv_metrics.rds"))
  