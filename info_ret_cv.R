library(here)

source(paste0(here(), "/load_data.R"))

# melting data
mdata <- load_mdata()

# interactors
interactors <- read_rds(paste0(here(), "/interactors.rds"))

# interactions
ints <- read_rds(paste0(here(), "/cv_ints.rds"))

# metrics trained on cross-validation folds of interactions
cv_metrics <- read_rds(paste0(here(), "/cv_metrics.rds"))
ITML <- cv_metrics$ITML
MMC_diag_neg <- cv_metrics$MMC_diag_neg

# Make test collection
make_test_collection <- function(q, rel_doc, num_docs) {
  
  q_interactors <- interactors %>% filter(Protein == q) %>% pull(Interactors) %>% unlist(.)
  
  docs <- c(rel_doc)
  ndocs <- 1
  while(ndocs < num_docs) { 
    
    doc <- sample(mdata$Protein, 1)
    if ( !(doc %in% c(q_interactors, docs)) ) {
      
      docs <- c(docs, doc)
      ndocs <- ndocs + 1
    }
  }
  
  return(list(Q = q, D = docs, R = tibble(Query = q, Document = rel_doc)))
}

IR_Mahal <- function(tc, M = NULL) {
  
  if (is.null(M))
    return(NA)
  
  # Get melting data of query
  q <- tibble(Protein = tc$Q) %>% 
    inner_join(mdata, by = "Protein") %>%
    select(T40:T64) %>%
    unlist(., use.names = F)
  
  # Compute the Mahalanobis distance between the query and every document.
  # Rank the documents in order of increasing Euclidean distance.
  dists <- tibble(Protein = tc$D) %>%
    inner_join(mdata, by = "Protein") %>%
    select(Protein, T40:T64) %>%
    nest(Sols = -Protein) %>%
    mutate(Sols = map(Sols, ~ unlist(.x, use.names = F))) %>%
    # Mahal. dist. squared
    mutate(Dist = map_dbl(Sols, ~ mahalanobis(x = q, center = .x, cov = M, inverted = T))) %>%
    arrange(Dist)
  
  # Extract the proteins
  rk <- dists$Protein
  
  # Ranks of relevant document
  return(match(tc$R$Document, rk))
}

# IR-Pear: Rank documents in order of increasing Pearson dissimilarity
# Arguments:
#   1. tc: test collection
IR_Pear <- function(tc) {
  
  # Get melting data of query
  q <- tibble(Protein = tc$Q) %>% 
    inner_join(mdata, by = "Protein") %>%
    select(T40:T64) %>%
    unlist(., use.names = F)
  
  # Compute the Pearson dissimilarity between the query and every document.
  # Rank the documents in order of increasing Pearson dissimilarity.
  cors <- tibble(Protein = tc$D) %>%
    inner_join(mdata, by = "Protein") %>%
    select(Protein, T40:T64) %>%
    nest(Sols = -Protein) %>%
    mutate(Sols = map(Sols, ~ unlist(.x, use.names = F))) %>%
    mutate(Cor = map_dbl(Sols, ~ (1 - cor(.x, q, method = "pearson")) / 2)) %>%
    arrange(Cor)
  
  # Extract the proteins
  rk <- cors$Protein
  
  # Ranks of relevant document
  return(match(tc$R$Document, rk))
}

set.seed(27)

# Test collections where ProteinA is the query and ProteinB is the relevant document
pb <- progress_bar$new(total = nrow(ints))
ints <- ints %>%
  mutate(TCA = map2(ProteinA, ProteinB, ~ {pb$tick(); make_test_collection(.x, .y, 20)}))

# Test collections where ProteinB is the query and ProteinA is the relevant document
pb <- progress_bar$new(total = nrow(ints))
ints <- ints %>%
  mutate(TCB = map2(ProteinB, ProteinA, ~ {pb$tick(); make_test_collection(.x, .y, 20)}))

# inds <- sample.int(nrow(ints), 100)

print("Running IR-ITML for each test collection...")
pb <- progress_bar$new(total = nrow(ints))
# pb <- progress_bar$new(total = 100)
res <- tibble(ITML_A = map2_int(ints$TCA, ints$Label, ~ {pb$tick(); IR_Mahal(.x, M = ITML[[.y]])}))

pb <- progress_bar$new(total = nrow(ints))
# pb <- progress_bar$new(total = 100)
res <- res %>%
  mutate(ITML_B = map2_int(ints$TCB, ints$Label, ~ {pb$tick(); IR_Mahal(.x, M = ITML[[.y]])}))

print("Running IR-MMC-diag-neg for each test collection...")
pb <- progress_bar$new(total = nrow(ints))
# pb <- progress_bar$new(total = 100)
res <- res %>%
  mutate(MMC_diag_neg_A = map2_int(ints$TCA, ints$Label, ~ {pb$tick(); IR_Mahal(.x, M = MMC_diag_neg[[.y]])}))

pb <- progress_bar$new(total = nrow(ints))
# pb <- progress_bar$new(total = 100)
res <- res %>%
  mutate(MMC_diag_neg_B = map2_int(ints$TCB, ints$Label, ~ {pb$tick(); IR_Mahal(.x, M = MMC_diag_neg[[.y]])}))

print("Running IR-Eucl for each test collection")
pb <- progress_bar$new(total = nrow(ints))
# pb <- progress_bar$new(total = 100)
res <- res %>%
  mutate(Eucl_A = map_int(ints$TCA, ~ {pb$tick(); IR_Mahal(.x, M = diag(9))}))

pb <- progress_bar$new(total = nrow(ints))
# pb <- progress_bar$new(total = 100)
res <- res %>%
  mutate(Eucl_B = map_int(ints$TCB, ~ {pb$tick(); IR_Mahal(.x, M = diag(9))}))

print("Running IR-Pear for each test collection")
pb <- progress_bar$new(total = nrow(ints))
# pb <- progress_bar$new(total = 100)
res <- res %>%
  mutate(Pear_A = map_int(ints$TCA, ~ {pb$tick(); IR_Pear(.x)}))

pb <- progress_bar$new(total = nrow(ints))
# pb <- progress_bar$new(total = 100)
res <- res %>%
  mutate(Pear_B = map_int(ints$TCB, ~ {pb$tick(); IR_Pear(.x)}))

write_rds(res, paste0(here(), "/IR_cv.rds"))

