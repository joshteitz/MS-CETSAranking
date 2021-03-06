library(here)

source(paste0(here(), "/load_data.R"))

# melting data
mdata <- load_mdata()
# mparams <- param_mdata(mdata)

# test collections
test_collections <- read_rds(paste0(here(), "/test_collections.rds"))

# ITML metrics trained on positive pairs
itml_pos <- read_rds(paste0(here(), "/itml_pos.rds"))

# Load test collections and learned metrics
# training_sets <- read_rds(paste0(here(), "/training_sets.rds"))
# learned_metrics <- read_rds(paste0(here(), "/learned_metrics.rds"))
# itml_metrics <- read_rds(paste0(here(), "/itml_metrics.rds"))
# mmc_metrics <- read_rds(paste0(here(), "/mmc_metrics.rds"))

# IR-Mahal: Rank documents in order of increasing Mahalanobis distance
# Arguments:
#   1. tc: test collection
#   2. M: matrix, default is diagonal matrix (Euclidean distance)
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
  return(match(tc$R$D, rk))
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
  return(match(tc$R$D, rk))
}

# IR-Par: Rank documents in order of increasing Pearson dissimilarity
# Arguments:
#   1. tc: test collection
IR_Par <- function(tc, M = NULL) {
  
  if (is.null(M))
    return(NA)
  
  # Get melting data of query and documents, and then standardize
  # each parameter value.
  q_docs <- bind_rows(
    tibble(Protein = tc$Q) %>% inner_join(mparams, by = "Protein"),
    tibble(Protein = tc$D) %>% inner_join(mparams, by = "Protein")
  ) %>% 
    make_unit_var(numeric_cols = c("Param_b", "Param_c", "Param_e"))
  
  # query with parametric melting data.
  q <- q_docs[1,] %>% select(starts_with("Param")) %>% unlist(., use.names = F)
  
  # Compute the Euclidean distance between the query and every document.
  # Rank the documents in order of increasing Parametric Euclidean distance.
  dists <- q_docs[-1,] %>%
    select(Protein, starts_with("Param")) %>%
    nest(Params = -Protein) %>%
    mutate(Params = map(Params, ~ unlist(.x, use.names = F))) %>%
    # Mahal. dist. squared
    mutate(Dist = map_dbl(Params, ~ mahalanobis(x = q, center = .x, cov = M, inverted = T))) %>%
    arrange(Dist)
  
  # Extract the proteins
  rk <- dists$Protein
  
  # Ranks of relevant documents
  return(match(tc$R$D, rk))
}

print("Run IR-Eucl on each test collection...")
pb <- progress_bar$new(total = nrow(test_collections))
eucl_res <- test_collections %>%
  mutate(Rank = map_int(TC, ~ {pb$tick(); IR_Mahal(.x, M = diag(1L, 9))})) %>%
  select(-TC)

print("Run IR-ITML-pos on each test collection...")
pb <- progress_bar$new(total = nrow(test_collections))
itml_pos_res <- test_collections %>%
  inner_join(itml_pos, by = c("Query" = "Query", "Rel_doc" = "Rel_doc")) %>%
  mutate(Rank = map2_int(TC, Metric, ~ {pb$tick(); IR_Mahal(.x, M = .y)})) %>%
  select(-TC, -Metric)

write_rds(eucl_res, paste0(here(), "/eucl_res.rds"))
write_rds(itml_pos_res, paste0(here(), "/itml_pos_res.rds"))

# set.seed(27)
# inds <- sample.int(nrow(training_sets), 100)

# # IR1: Results for IR-Eucl, IR-Pear, and IR-Par
# IR1 <- tibble(Q = itml_metrics$Q, R_doc = itml_metrics$R_doc)
# 
# print("Results for IR-Eucl...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR1 <- IR1 %>% mutate(Eucl = map_int(training_sets$TC, ~{pb$tick(); IR_Mahal(.x, diag(9))}))
# 
# print("Results for IR-Pear...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR1 <- IR1 %>% mutate(Pear = map_int(training_sets$TC, ~{pb$tick(); IR_Pear(.x)}))
# 
# print("Results for IR-Par...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR1 <- IR1 %>% mutate(Par = map_int(training_sets$TC, ~{pb$tick(); IR_Par(.x, diag(3))}))
# 
# # IR2: Results for IR-ITML, IR-ITML-eg, IR-ITML-Par, IR-ITML-Par-neg
# IR2 <- tibble(Q = itml_metrics$Q, R_doc = itml_metrics$R_doc)
# 
# print("Results for IR-ITML...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR2 <- IR2 %>% mutate(ITML = map2_int(training_sets$TC, itml_metrics$ITML, ~{pb$tick(); IR_Mahal(.x, .y)}))
# 
# print("Results for IR-ITML-neg...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR2 <- IR2 %>% mutate(ITML_neg = map2_int(training_sets$TC, itml_metrics$ITML_neg, ~{pb$tick(); IR_Mahal(.x, .y)}))
# 
# print("Results for IR-ITML-Par...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR2 <- IR2 %>% mutate(ITML_par = map2_int(training_sets$TC, itml_metrics$ITML_par, ~{pb$tick(); IR_Par(.x, .y)}))
# 
# print("Results for IR-ITML-Par-neg...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR2 <- IR2 %>% mutate(ITML_par_neg = map2_int(training_sets$TC, itml_metrics$ITML_par_neg, ~{pb$tick(); IR_Par(.x, .y)}))

# # IR: Results for IR-Eucl, IR-Pear, IR-Par, IR-ITML, IR-ITML-neg, IR-MMC-neg
# IR <- tibble(Q = learned_metrics$Q, R_doc = learned_metrics$R_doc)
# 
# print("Results for IR-Eucl...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR <- IR %>% mutate(Eucl = map_int(training_sets$TC, ~{pb$tick(); IR_Mahal(.x, diag(9))}))
# 
# print("Results for IR-Pear...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR <- IR %>% mutate(Pear = map_int(training_sets$TC, ~{pb$tick(); IR_Pear(.x)}))
# 
# print("Results for IR-Par...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR <- IR %>% mutate(Par = map_int(training_sets$TC, ~{pb$tick(); IR_Par(.x, diag(3))}))
# 
# print("Results for IR-ITML...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR <- IR %>% mutate(ITML = map2_int(training_sets$TC, learned_metrics$ITML, ~{pb$tick(); IR_Mahal(.x, .y)}))
# 
# print("Results for IR-ITML-neg...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR <- IR %>% mutate(ITML_neg = map2_int(training_sets$TC, learned_metrics$ITML_neg, ~{pb$tick(); IR_Mahal(.x, .y)}))
# 
# print("Results for IR-MMC-id...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR <- IR %>% mutate(MMC_neg = map2_int(training_sets$TC, learned_metrics$MMC_neg, ~{pb$tick(); IR_Mahal(.x, .y)}))

# # IR3: Results for IR-MMC-id, IR-MMC-rand, IR-MMC-Par-id, IR-MMC-Par-rand
# IR3 <- tibble(Q = mmc_metrics$Q, R_doc = mmc_metrics$R_doc)
# 
# print("Results for IR-MMC-id...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR3 <- IR3 %>% mutate(MMC_id = map2_int(training_sets$TC, mmc_metrics$MMC_id, ~{pb$tick(); IR_Mahal(.x, .y)}))
# 
# print("Results for IR-MMC-rand...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR3 <- IR3 %>% mutate(MMC_rand = map2_int(training_sets$TC, mmc_metrics$MMC_rand, ~{pb$tick(); IR_Mahal(.x, .y)}))
# 
# print("Results for IR-MMC-Par-id...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR3 <- IR3 %>% mutate(MMC_par_id = map2_int(training_sets$TC, mmc_metrics$MMC_par_id, ~{pb$tick(); IR_Par(.x, .y)}))
# 
# print("Results for IR-MMC-Par-rand...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# IR3 <- IR3 %>% mutate(MMC_par_rand = map2_int(training_sets$TC, mmc_metrics$MMC_par_rand, ~{pb$tick(); IR_Par(.x, .y)}))

# write_rds(IR1, paste0(here(), "/IR1.rds"))
# write_rds(IR2, paste0(here(), "/IR2.rds"))
# write_rds(IR, paste0(here(), "/IR.rds"))

# IR1 <- read_rds(paste0(here(), "/IR1.rds"))
# IR2 <- read_rds(paste0(here(), "/IR2.rds"))
# 
# # Mean ranks
# IR1 %>% summarize(across(Eucl:Par, ~ mean(.x, na.rm = T)))
# IR2 %>% summarize(across(ITML:ITML_par_neg, ~ mean(.x, na.rm = T)))
# 
# # Indicator for test collections that result in a null metric
# null_metrics <- map_lgl(itml_metrics$ITML_par_neg, ~ is.null(.x))
# 
# IR1[!null_metrics,] %>% summarize(across(Eucl:Par, ~ mean(.x, na.rm = T)))
# 
# # of positive examples in each training set
# num_pos <- map_int(training_sets$TC, ~ length(.x$Pos))
# 
# # correlations
# rks <- IR1$Eucl
# cor_df <- tibble(
#   Num_pos = num_pos,
#   Rank = rks
# ) %>%
#   filter(Num_pos >= 5)
# 
# cor_df %>%
#   group_by(Rng = cut(Num_pos, breaks = seq(5, 120, by = 5), right = F)) %>%
#   summarize(Mean_rank = mean(Rank))
#   
# 
# # mean_ranks: find the mean rank for each IR system
# # Arguments:
# #   1. IR_res: IR system results
# #   2. IR_cols: specify IR systems
# mean_ranks <- function(IR_res, IR_cols = c("Eucl", "Pear", "Par")) {
#   return(IR_res %>% summarize(across({{IR_cols}}, mean)))
# }
