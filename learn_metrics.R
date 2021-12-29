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

set.seed(27)
inds <- sample.int(nrow(training_sets), 100)

print("Learning ITML metrics...")
pb <- progress_bar$new(total = nrow(training_sets))
pb <- progress_bar$new(total = 100)
itml_metrics <- tibble(ITML = map(trs[inds], ~ {pb$tick(); learn_itml(.x, min_num = 10, use_neg_pairs = F)}))

print("Learning ITML metrics with negative training examples...")
pb <- progress_bar$new(total = nrow(training_sets))
pb <- progress_bar$new(total = 100)
itml_metrics <- itml_metrics %>%
  mutate(ITML_neg = map(trs[inds], ~ {pb$tick(); learn_itml(.x, min_num = 10, use_neg_pairs = T)}))

print("Learning ITML metrics on parametric training data")
pb <- progress_bar$new(total = nrow(training_sets))
pb <- progress_bar$new(total = 100)
itml_metrics <- itml_metrics %>%
  mutate(ITML_par = map(trs_p[inds], ~ {pb$tick(); learn_itml(.x, min_num = 10, use_neg_pairs = F)}))

print("Learning ITML metrics on parametric training data with negative training examples")
pb <- progress_bar$new(total = nrow(training_sets))
pb <- progress_bar$new(total = 100)
itml_metrics <- itml_metrics %>%
  mutate(ITML_par_neg = map(trs_p[inds], ~ {pb$tick(); learn_itml(.x, min_num = 10, use_neg_pairs = T)}))

write_rds(itml_metrics, paste0(here(), "/itml_metrics.rds"))

# i <- sample.int(nrow(training_sets), 1) # 22913
# tr <- training_sets$Tr[[i]]
# tr
# r_to_py(tr)
# 
# learn_mmc(tr, min_num = 10, use_neg_pairs = F)

