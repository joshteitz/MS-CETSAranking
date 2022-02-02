library(reticulate)
library(tidyverse)
library(conflicted)
library(progress)
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

# use_condaenv("metric-learning")

source(paste0(here(), "/load_data.R"))

source_python("learn_metrics.py")

# load test collections
test_collections <- read_rds(paste0(here(), "/test_collections.rds"))

# load training data
tr_pos <- read_rds(paste0(here(), "/tr_pos.rds"))

# Function to remove any training pair containing a particular document
rm_tr_pair <- function(tr, doc) {
  # find row number of `doc`
  rnum <- which(rownames(tr$X) == doc)
  # remove from `pairs_indices` and `y_pairs`
  pairs_indices <- tr$pairs_indices[-(rnum-1)]
  y_pairs <- tr$y_pairs[-(rnum-1)]
  return(list(X = tr$X, pairs_indices = pairs_indices, y_pairs = y_pairs))
}

# For each test collection, learn an ITML metric based on the query's positive training pairs.
print("Learning an ITML metric based on each protein's positive training pairs...")
pb <- progress_bar$new(total = nrow(test_collections))
itml_pos <- test_collections %>%
  inner_join(tr_pos, by = c("Query" = "Protein")) %>%
  mutate(Metric = map2(Rel_doc, Tr, ~ {
    pb$tick();
    tr <- rm_tr_pair(.y, .x);
    learn_itml(tr)
  })) %>%
  select(-TC, -Tr)

write_rds(itml_pos, paste0(here(), "/itml_pos.rds"))





# # import training sets
# training_sets <- read_rds(paste0(here(), "/training_sets.rds"))
# 
# # training sets
# trs <- training_sets$Tr
# 
# # set.seed(27)
# # inds <- sample.int(nrow(training_sets), 100)
#  
# # ITML metrics
# learned_metrics <- tibble(
#   Q = map_chr(training_sets$TC, ~ .x$Q),
#   R_doc = map_chr(training_sets$TC, ~ .x$R$Document),
#   Tr = trs,
# ) %>%
#   mutate(Tr_size = map_int(Tr, ~ .x$y_pairs %>% length))
# 
# # itml_metrics <- itml_metrics[inds,]
# 
# # Set minimum number of training examples
# MIN_NUM = 6
# 
# print("Learning ITML metrics...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# learned_metrics <- learned_metrics %>%
#   mutate(ITML = pmap(list(Q, R_doc, Tr), ~ {pb$tick(); learn_itml(..3, min_num = MIN_NUM, use_neg_pairs = F)}))
# 
# print("Learning ITML metrics with negative training examples...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# learned_metrics <- learned_metrics %>%
#   mutate(ITML_neg = pmap(list(Q, R_doc, Tr), ~ {pb$tick(); learn_itml(..3, min_num = MIN_NUM, use_neg_pairs = T)}))
# 
# print("Learning ITML metrics on parametric training data...")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# itml_metrics <- itml_metrics %>%
#   mutate(ITML_par = pmap(list(Q, R_doc, Tr_p), ~ {pb$tick(); learn_itml(..3, min_num = MIN_NUM, use_neg_pairs = F)}))
# 
# print("Learning ITML metrics on parametric training data with negative training examples")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# itml_metrics <- itml_metrics %>%
#   mutate(ITML_par_neg = pmap(list(Q, R_doc, Tr_p), ~ {pb$tick(); learn_itml(..3, min_num = MIN_NUM, use_neg_pairs = T)}))

# # MMC metrics
# mmc_metrics <- tibble(
#   Q = map_chr(training_sets$TC, ~ .x$Q),
#   R_doc = map_chr(training_sets$TC, ~ .x$R$Document),
#   Tr = trs,
#   Tr_p = trs_p
# )

# print("Learning MMC metrics with identity initialization and negative training examples")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# learned_metrics <- learned_metrics %>%
#   mutate(MMC_neg = pmap(list(Q, R_doc, Tr), ~ {
#     pb$tick(); 
#     learn_mmc(..3, min_num = MIN_NUM, use_neg_pairs = T, diag = T, initialization = "identity")
#   }))

# print("Learning MMC metrics with random initialization and negative training examples")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# mmc_metrics <- mmc_metrics %>%
#   mutate(MMC_rand = pmap(list(Q, R_doc, Tr), ~ {
#     pb$tick(); 
#     learn_mmc(..3, min_num = MIN_NUM, use_neg_pairs = T, diag = T, initialization = "random")
#   }))

# print("Learning MMC metrics on parametric data with identity initialization and negative training examples")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# mmc_metrics <- mmc_metrics %>%
#   mutate(MMC_par_id = pmap(list(Q, R_doc, Tr_p), ~ {
#     pb$tick(); 
#     learn_mmc(..3, min_num = MIN_NUM, use_neg_pairs = T, diag = T, initialization = "identity")
#   }))
# 
# print("Learning MMC metrics on parametric data with random initialization and negative training examples")
# pb <- progress_bar$new(total = nrow(training_sets))
# # pb <- progress_bar$new(total = 100)
# mmc_metrics <- mmc_metrics %>%
#   mutate(MMC_par_rand = pmap(list(Q, R_doc, Tr_p), ~ {
#     pb$tick(); 
#     learn_mmc(..3, min_num = MIN_NUM, use_neg_pairs = T, diag = T, initialization = "random")
#   }))

# write_rds(itml_metrics, paste0(here(), "/itml_metrics.rds"))

# write_rds(learned_metrics, paste0(here(), "/learned_metrics.rds"))

# set.seed(27)
# i <- sample.int(nrow(training_sets), 1) # 22913
# tr <- training_sets$Tr[[i]]
# tr
# r_to_py(tr)
# 
# res <- learn_mmc(tr, min_num = 10, use_neg_pairs = T, diag = T, initialization = "identity")

# tr = (itml_metrics %>% filter(Q == "P51668" & R_doc == "O15033"))$Tr[[1]]
# res <- learn_itml(tr, min_num = 10, use_neg_pairs = T)
# 
# py_save_object(tr, "tr.pickle")
