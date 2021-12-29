library(here)

source(paste0(here(), "/load_data.R"))

 # melting data from Table S7
mdata <- load_mdata()

# # convert melting data to parametric melting data
mparams <- param_mdata(mdata) %>% drop_na()

# load test collections
test_collections <- read_rds(paste0(here(), "/test_collections.rds"))

# Number of relevant documents for each query
test_collections <- test_collections %>% mutate(Num_rel = map_int(Rel_documents, length))

# make_training_set: make training set out of test collection.
# Arguments:
#   1. tc: test collecction
#   2. parametric: should training set be based on parametric melting data
make_training_set <- function(tc, parametric = F) {
  
  # if all relevant documents are part of the test collection
  if ( length(tc$Pos) == 0 )
    return(NULL)
  
  # Proteins
  X <- tibble(Protein = c(tc$Q, tc$Pos, tc$Neg))
  
  # Melting data for proteins
  if (parametric) {
    X <- X %>% inner_join(mparams, by = "Protein") %>%
      make_unit_var(numeric_cols = c(Param_b, Param_c, Param_e))
  } else {
    X <- X %>% inner_join(mdata, by = "Protein") %>% select(-T37)
  }
  
  # Convert to matrix
  X <- X %>% select(-Protein) %>% as.matrix()
  colnames(X) <- NULL
  
  # Convert to list of vectors
  X <- array_tree(X, margin = 1)
  
  # positive pairs and negative pairs
  pairs_indices <- expand.grid(0, 1:(length(X) - 1)) %>% as.matrix()
  colnames(pairs_indices) <- NULL
  pairs_indices <- array_tree(pairs_indices, margin = 1)
  pairs_indices <- map(pairs_indices, as.integer)
  
  # labels for pairs
  y_pairs <- c(rep(1L, length(tc$Pos)), rep(-1L, length(tc$Neg)))
  
  return(list(X = X, pairs_indices = pairs_indices, y_pairs = y_pairs))
}

# set.seed(27)
# # test collection with one relevant document
# NUM_REL = 1
# query <- (test_collections %>% filter(Num_rel == NUM_REL))$Query %>% sample(size = 1)
# tc <- (learned_metrics %>% filter(Query == query) %>% slice_sample(n = 1))$TC[[1]]
# 
# # test collection with two relevant documents
# NUM_REL = 2
# query <- (test_collections %>% filter(Num_rel == NUM_REL))$Query %>% sample(size = 1)
# tc1 <- (learned_metrics %>% filter(Query == query) %>% slice_sample(n = 1))$TC[[1]]
# 
# # test collection with three relevant documents
# NUM_REL = 3
# query <- (test_collections %>% filter(Num_rel == NUM_REL))$Query %>% sample(size = 1)
# tc2 <- (learned_metrics %>% filter(Query == query) %>% slice_sample(n = 1))$TC[[1]]
# 
# # test collection with four relevant documents
# NUM_REL = 4
# query <- (test_collections %>% filter(Num_rel == NUM_REL))$Query %>% sample(size = 1)
# tc3 <- (learned_metrics %>% filter(Query == query) %>% slice_sample(n = 1))$TC[[1]]
# 
# # test collection with five relevant documents
# NUM_REL = 5
# query <- (test_collections %>% filter(Num_rel == NUM_REL))$Query %>% sample(size = 1)
# tc4 <- (learned_metrics %>% filter(Query == query) %>% slice_sample(n = 1))$TC[[1]]

# Make a training set for each test collection
training_sets <- map2_dfr(test_collections$Query, test_collections$TCs, ~ tibble(Query = .x, TC = .y))
pb <- progress_bar$new(total = nrow(training_sets))
pb1 <- progress_bar$new(total = nrow(training_sets))
training_sets <- training_sets %>%
  mutate(Tr = map(TC, ~ {pb$tick(); make_training_set(.x)})) %>%
  mutate(Tr_param = map(TC, ~ {pb1$tick(); make_training_set(.x, parametric = T)}))

# # Size of each training set
# training_sets <- training_sets %>%
#   mutate(Tr_size = map_int(Training, ~ if_else(is.null(.x), 0L, length(.x$y_pairs))))

write_rds(training_sets, paste0(here(), "/training_sets.rds"))

