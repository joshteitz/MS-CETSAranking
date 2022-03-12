library(here)

source(here("load_data.R"))

SZ_NULL = 10000

# FUNCTION: IR-Eucl
IR_Eucl <- function(tc, mdata) {
  
  # melting data of query
  md_Q <- tibble(Protein = tc$Q) %>% 
    inner_join(mdata %>% select(-T37), by = c("Protein")) %>%
    nest(Sols = starts_with("T")) %>%
    mutate(Sols = map(Sols, ~ unlist(.x, use.names = F))) %>%
    pull(Sols) %>%
    unlist()
  
  # Return NULL in the case when query does not have melting data
  if (is.null(md_Q)) return(NULL)
  
  # melting data of documents
  md_D <- tibble(Protein = tc$D) %>%
    inner_join(mdata %>% select(-T37), by = c("Protein")) %>%
    nest(Sols = starts_with("T")) %>%
    mutate(Sols = map(Sols, ~ unlist(.x, use.names = F)))
  
  # Return Null in the case when documents have no melting data
  if (nrow(md_D) == 0) return(NULL)
  
  # Convert melting data of documents to matrix
  prots <- md_D$Protein
  md_D <- matrix(unlist(md_D$Sols), byrow = T, nrow = length(md_D$Sols))
  
  # Compute (squared) Euclidean distance between query and each document
  dists <- mahalanobis(md_D, md_Q, cov = diag(1, length(md_Q)), inverted = T)
  names(dists) <- prots
  
  # Rank the proteins in order of increasing Euclidean distance from the query.
  rk <- dists %>% sort %>% names
  
  # return the ranks of the relevant documents
  rel_docs <- tc$R$D
  
  # Ranks of relevant documents
  rk_rel <- which(rk %in% rel_docs)
  
  if (length(rk_rel) == 0) {
    return(NULL)
  } else {
    return(rk_rel)
  }
}

# FUNCTION: IR-Par
IR_Par <- function(tc, mparams) {
  
  # melting data of query
  md_Q <- tibble(Protein = tc$Q) %>% inner_join(mparams, by = "Protein")
  
  # Return NULL in the case when query does not have melting data
  if (nrow(md_Q) == 0) return(NULL)
  
  # melting data of documents
  md_D <- tibble(Protein = tc$D ) %>% inner_join(mparams, by = "Protein")
  
  # Return NULL in the case when documents have no melting data
  if (is.null(md_Q)) return(NULL)
  
  # Combine and standardize melting data of query and documents
  md <- bind_rows(md_Q, md_D) %>%
    make_unit_var(., numeric_cols = c("Param_b", "Param_c", "Param_e")) %>%
    nest(Sols = starts_with("Param")) %>%
    mutate(Sols = map(Sols, ~ unlist(.x, use.names = F)))
  
  # melting data of query
  md_Q <- md[1,]$Sols[[1]]
  
  # melting data of docuements
  md_D <- matrix(unlist(md[-1,]$Sols), byrow = T, nrow = length(md[-1,]$Sols))
  prots <- md[-1,]$Protein
  
  # Compute (squared) Euclidean distance between query and each document
  dists <- mahalanobis(md_D, md_Q, cov = diag(1, length(md_Q)), inverted = T)
  names(dists) <- prots
  
  # Rank the proteins in order of increasing Euclidean distance from the query.
  rk <- dists %>% sort %>% names
  
  # return the ranks of the relevant documents
  rel_docs <- tc$R$D
  
  # Ranks of relevant documents
  rk_rel <- which(rk %in% rel_docs)
  
  if (length(rk_rel) == 0) {
    return(NULL)
  } else {
    return(rk_rel)
  }
}

# FUNCTION: compute the Pearson dissimilarity between two vectors
pear_dissim <- function(x, y) {
  (1 - cor(x, y, method = "pearson")) / 2
}

# FUNCTION: IR-Pear
IR_Pear <- function(tc, mdata) {
  
  # melting data of query
  md_Q <- tibble(Protein = tc$Q) %>% 
    inner_join(mdata %>% select(-T37), by = c("Protein")) %>%
    nest(Sols = starts_with("T")) %>%
    mutate(Sols = map(Sols, ~ unlist(.x, use.names = F))) %>%
    pull(Sols) %>%
    unlist()
  
  # Return NULL in the case when query does not have melting data
  if (is.null(md_Q)) return(NULL)
  
  # melting data of documents
  md_D <- tibble(Protein = tc$D) %>%
    inner_join(mdata %>% select(-T37), by = c("Protein")) %>%
    nest(Sols = starts_with("T")) %>%
    mutate(Sols = map(Sols, ~ unlist(.x, use.names = F)))
  
  # Return Null in the case when documents have no melting data
  if (nrow(md_D) == 0) return(NULL)
  
  # Compute the Pearson dissimilarity between query and each document
  dists <- map_dbl(md_D$Sols, ~ pear_dissim(x = .x, y = md_Q))
  names(dists) <- md_D$Protein
  
  # Rank the proteins in order of increasing Euclidean distance from the query.
  rk <- dists %>% sort %>% names
  
  # return the ranks of the relevant documents
  rel_docs <- tc$R$D
  
  # Ranks of relevant documents
  rk_rel <- which(rk %in% rel_docs)
  
  if (length(rk_rel) == 0) {
    return(NULL)
  } else {
    return(rk_rel)
  }
}

# FUNCTION: compute average precision based on the ranks of relevant documents
#   rks: ranks of relevant documents in order of increasing rank
avg_prec <- function(rks) {
  num_rel <- length(rks)
  ap <- sum(1:num_rel / rks) / num_rel
  return(ap)
}

# FUNCTION: Compute an empirical null distribution of average precision
#  num_rel: number of relevant documents
#  num: total number of documents
#  sz_null: the size of the empirical null distribution
null_avg_prec <- function(num_rel, num, sz_null) {
  
  vals <- vector(mode = "double", length = sz_null)
  for (i in 1:sz_null) {
    ranks <- sample.int(num, num_rel) %>% sort
    vals[i] <- avg_prec(ranks)
  }
  
  return(vals)
}

# FUNCTION: Perform IR-Eucl, IR-Par, IR-Pear on test collections.
#   tcs: test collections
#   sz_null: size of empirical null distribution
perform_IR <- function(tcs, mdata, mparams, sz_null) {
  
  print("Running IR-Eucl on test collections...")
  pb <- progress_bar$new(total = nrow(tcs))
  ir_eucl <- tcs %>% 
    mutate(Ranking = map(TC, ~ {pb$tick(); IR_Eucl(.x, mdata)} )) %>%
    select(-TC) %>%
    drop_na(Ranking)
  
  print("Computing average precision scores for IR-Eucl rankings...")
  pb <- progress_bar$new(total = nrow(ir_eucl))
  ir_eucl <- ir_eucl %>% mutate(Avg_prec = map_dbl(Ranking, ~ {pb$tick(); avg_prec(.x)}))
  
  print("Running IR-Pear on test collections...")
  pb <- progress_bar$new(total = nrow(tcs))
  ir_pear <- tcs %>% 
    mutate(Ranking = map(TC, ~ {pb$tick(); IR_Pear(.x, mdata)} )) %>%
    select(-TC) %>%
    drop_na(Ranking)
  
  print("Computing average precision scores for IR-Pear rankings...")
  pb <- progress_bar$new(total = nrow(ir_pear))
  ir_pear <- ir_pear %>% mutate(Avg_prec = map_dbl(Ranking, ~ {pb$tick(); avg_prec(.x)}))
  
  print("Running IR-Par on test collections...")
  pb <- progress_bar$new(total = nrow(tcs))
  ir_par <- tcs %>% 
    mutate(Ranking = map(TC, ~ {pb$tick(); IR_Par(.x, mparams)})) %>%
    select(-TC) %>%
    drop_na(Ranking)
  
  print("Computing average precision scores for IR-Par rankings...")
  pb <- progress_bar$new(total = nrow(ir_par))
  ir_par <- ir_par %>% mutate(Avg_prec = map_dbl(Ranking, ~ {pb$tick(); avg_prec(.x)}))
  
  print(
    "Computing empirical null distribution of average precision for each test collection..."
  )
  pb <- progress_bar$new(total = nrow(tcs))
  tcs <- tcs %>% mutate(Null_distr = map(TC, ~ {
    pb$tick();
    null_avg_prec(num_rel = .x$R %>% nrow, num = .x$D %>% length, sz_null)
  }))
  
  print("Computing an empirical p-value for each IR-Eucl average precision score...")
  ir_eucl <- ir_eucl %>%
    inner_join(tcs %>% select(-TC), by = "Bait") %>%
    mutate(Pval = map2_dbl(Avg_prec, Null_distr, ~ sum(.y >= .x) / length(.y))) %>%
    select(-Null_distr)
  
  print("Computing an empirical p-value for each IR-Pear average precision score...")
  ir_pear <- ir_pear %>%
    inner_join(tcs %>% select(-TC), by = "Bait") %>%
    mutate(Pval = map2_dbl(Avg_prec, Null_distr, ~ sum(.y >= .x) / length(.y))) %>%
    select(-Null_distr)
  
  print("Computing an empirical p-value for each IR-Par average precision score...")
  ir_par <- ir_par %>%
    inner_join(tcs %>% select(-TC), by = "Bait") %>%
    mutate(Pval = map2_dbl(Avg_prec, Null_distr, ~ sum(.y >= .x) / length(.y))) %>%
    select(-Null_distr)
  
  # null distribution average precision scores in matrix form
  ap_scores <- tcs$Null_distr
  ap_scores <- matrix(unlist(ap_scores), byrow = T, nrow = length(ap_scores))
  
  # Null distribution of mAP
  null_map <- colMeans(ap_scores)
  
  # mAP values
  map_eucl <- ir_eucl$Avg_prec %>% mean
  map_pear <- ir_pear$Avg_prec %>% mean
  map_par <- ir_par$Avg_prec %>% mean
  
  # p-vals
  pval_eucl <- sum(null_map >= map_eucl) / length(null_map)
  pval_pear <- sum(null_map >= map_pear) / length(null_map)
  pval_par <- sum(null_map >= map_par) / length(null_map)
  
  res <- list(
    IR_Eucl = list(rankings = ir_eucl, mAP = map_eucl, pval = pval_eucl),
    IR_Pear = list(rankings = ir_pear, mAP = map_pear, pval = pval_pear),
    IR_Par = list(rankings = ir_par, mAP = map_par, pval = pval_par)
  )
  
  return(res)
}

# Load melting data from three cell lines
mdata_k562 <- load_mdata("K562")
mdata_293t <- load_mdata("HEK293T")
mdata_hct116 <- load_mdata("HCT116")

# Make parametric melting data
mparams_k562 <- param_mdata(mdata_k562)
mparams_293t <- param_mdata(mdata_293t)
mparams_hct116 <- param_mdata(mdata_hct116)

# Melting data: K562
# Bait/prey data: 293T
# Interactions: S2 with more than two citations
tcs1 <- read_rds(here("tcs", "tcs1.rds"))
ir1 <- perform_IR(tcs1, mdata_k562, mparams_k562, SZ_NULL)
write_rds(ir1, here("ir", "ir1.rds"))

# Melting data: K562
# Bait/prey data: 293T
# Interactions: 293T
tcs2 <- read_rds(here("tcs", "tcs2.rds"))
ir2 <- perform_IR(tcs2, mdata_k562, mparams_k562, SZ_NULL)
write_rds(ir2, here("ir", "ir2.rds"))

# Melting data: K562
# Bait/prey data: HCT116
# Interactions: S2 with more than two citations
tcs3 <- read_rds(here("tcs", "tcs3.rds"))
ir3 <- perform_IR(tcs3, mdata_k562, mparams_k562, SZ_NULL)
write_rds(ir3, here("ir", "ir3.rds"))

# Melting data: K562
# Bait/prey data: HCT116
# Interactions: HCT116
tcs4 <- read_rds(here("tcs", "tcs4.rds"))
ir4 <- perform_IR(tcs4, mdata_k562, mparams_k562, SZ_NULL)
write_rds(ir4, here("ir", "ir4.rds"))

# Melting data: 293T
# Bait/prey data: 293T
# Interactions: S2
tcs1 <- read_rds(here("tcs", "tcs1.rds"))
ir5 <- perform_IR(tcs1, mdata_293t, mparams_293t, SZ_NULL)
write_rds(ir5, here("ir", "ir5.rds"))

# Melting data: 293T
# Bait/prey data: 293T
# Interactions: 293T
tcs2 <- read_rds(here("tcs", "tcs2.rds"))
ir6 <- perform_IR(tcs2, mdata_293t, mparams_293t, SZ_NULL)
write_rds(ir6, here("ir", "ir6.rds"))

# Melting data: HCT116
# Bait/prey data: HCT116
# Interactions: S2
tcs3 <- read_rds(here("tcs", "tcs3.rds"))
ir7 <- perform_IR(tcs3, mdata_hct116, mparams_hct116, SZ_NULL)
write_rds(ir7, here("ir", "ir7.rds"))

# Melting data: HCT116
# Bait/prey data: HCT116
# Interactions: HCT116
tcs4 <- read_rds(here("tcs", "tcs4.rds"))
ir8 <- perform_IR(tcs4, mdata_hct116, mparams_hct116, SZ_NULL)
write_rds(ir8, here("ir", "ir8.rds"))