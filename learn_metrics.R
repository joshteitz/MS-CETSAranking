library(here)
library(reticulate)
library(tidyverse)
library(conflicted)
library(progress)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

# use_condaenv("metric-learning")

source_python("learn_metrics.py")

# import training sets
training_sets <- read_rds(paste0(here(), "/training_sets.rds"))

# training sets
trs <- training_sets$Tr
trs_p <- training_sets$Tr_param

# set.seed(27)
# inds <- sample.int(nrow(training_sets), 100)

itml_metrics <- tibble(
  Q = map_chr(training_sets$TC, ~ .x$Q),
  R_doc = map_chr(training_sets$TC, ~ .x$R$Document),
  Tr = trs,
  Tr_p = trs_p
)

# itml_metrics <- itml_metrics[inds,]

print("Learning ITML metrics...")
pb <- progress_bar$new(total = nrow(training_sets))
# pb <- progress_bar$new(total = 100)
itml_metrics <- itml_metrics %>%
  mutate(ITML = pmap(list(Q, R_doc, Tr), ~ {pb$tick(); learn_itml(..3, min_num = 10, use_neg_pairs = F)}))

print("Learning ITML metrics with negative training examples...")
pb <- progress_bar$new(total = nrow(training_sets))
# pb <- progress_bar$new(total = 100)
itml_metrics <- itml_metrics %>%
  mutate(ITML_neg = pmap(list(Q, R_doc, Tr), ~ {pb$tick(); print(paste(..1, ..2)); learn_itml(..3, min_num = 10, use_neg_pairs = T)}))

print("Learning ITML metrics on parametric training data...")
pb <- progress_bar$new(total = nrow(training_sets))
# pb <- progress_bar$new(total = 100)
itml_metrics <- itml_metrics %>%
  mutate(ITML_par = pmap(list(Q, R_doc, Tr_p), ~ {pb$tick(); learn_itml(..3, min_num = 10, use_neg_pairs = F)}))

print("Learning ITML metrics on parametric training data with negative training examples")
pb <- progress_bar$new(total = nrow(training_sets))
# pb <- progress_bar$new(total = 100)
itml_metrics <- itml_metrics %>%
  mutate(ITML_par_neg = pmap(list(Q, R_doc, Tr_p), ~ {pb$tick(); learn_itml(..3, min_num = 10, use_neg_pairs = T)}))

write_rds(itml_metrics, paste0(here(), "/itml_metrics.rds"))

# i <- sample.int(nrow(training_sets), 1) # 22913
# tr <- training_sets$Tr[[i]]
# tr
# r_to_py(tr)
# 
# learn_mmc(tr, min_num = 10, use_neg_pairs = F)

# tr = (itml_metrics %>% filter(Q == "P51668" & R_doc == "O15033"))$Tr[[1]]
# res <- learn_itml(tr, min_num = 10, use_neg_pairs = T)
# 
# py_save_object(tr, "tr.pickle")
