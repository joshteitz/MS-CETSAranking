library(here)

source(paste0(here(), "/load_data.R"))

# melting data from Table S7
mdata <- load_mdata()

# convert melting data to parametric melting data
print("Converting melting data to parametric melting data...")
mparams <- param_mdata(mdata) %>% drop_na()

# pairwise interactions from table S2
ints_s2 <- load_interactions() %>% mutate(Num_pubs = Num_pubs %>% as.integer())

# Remove each interaction where at least one protein does not have parametric melting data
ints_s2 <- ints_s2 %>%
  semi_join(mparams, by = c("ProteinA" = "Protein")) %>%
  semi_join(mparams, by = c("ProteinB" = "Protein"))

# Ensure there are no duplicate interactions
ints_s2 <- ints_s2 %>%
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
ints_complex <- cdata %>%
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
ints_complex <- ints_complex %>%
  mutate(Proteins = map2(ProteinA, ProteinB, ~ sort(c(.x, .y)))) %>%
  mutate(ProteinA = map_chr(Proteins, ~ .x[1])) %>%
  mutate(ProteinB = map_chr(Proteins, ~ .x[2])) %>%
  distinct(ProteinA, ProteinB, .keep_all = T) %>%
  select(ProteinA, ProteinB, Complex)

# union of ints_s2 and ints_complex
ints_all <- bind_rows(ints_s2 %>% select(-Num_pubs), ints_complex %>% select(-Complex)) %>%
  distinct(ProteinA, ProteinB)

# Each protein's interactors based on interactions of ints_all
interactors_all <- bind_rows(
  ints_all,
  ints_all %>% select(ProteinA = ProteinB, ProteinB = ProteinA)
) %>%
  nest(Interactors = -ProteinA) %>%
  mutate(Interactors = map(Interactors, ~ .x$ProteinB %>% unique)) %>%
  select(Protein = ProteinA, Interactors)

# Each protein's interactors based on interactions with two or more publications
interactors_2pubs <- bind_rows(
  ints_s2 %>% filter(Num_pubs >= 2),
  ints_s2 %>% filter(Num_pubs >= 2) %>% select(ProteinA = ProteinB, ProteinB = ProteinA)
) %>%
  nest(Interactors = -ProteinA) %>%
  mutate(Interactors = map(Interactors, ~ .x$ProteinB %>% unique)) %>%
  select(Protein = ProteinA, Interactors)

# Random sample of `n` objects from vector `x`. The objects in `excl` are excluded from the random sample.
random_sample <- function(x, n, excl) {
  valid_objs <- x[!x %in% excl]
  return(sample(valid_objs, size = n))
}

# Make a test collection for each interaction with two or more publications.
test_collections <- bind_rows(
  ints_s2 %>% filter(Num_pubs >= 2) %>% select(Query = ProteinA, Rel_doc = ProteinB),
  ints_s2 %>% filter(Num_pubs >= 2) %>% select(Query = ProteinB, Rel_doc = ProteinA)
)
# Each test collection contains a single query, single relevant document and nineteen non-relevant documents.
NUM_NON_REL = 19
set.seed(27)
print("Making test collections...")
pb <- progress_bar$new(total = nrow(test_collections))
test_collections <- test_collections %>%
  inner_join(interactors_all, by = c("Query" = "Protein")) %>%
  mutate(TC = pmap(list(Query, Rel_doc, Interactors), ~{pb$tick(); list(
    Q = ..1,
    D = c(..2, random_sample(mparams$Protein, NUM_NON_REL, excl = ..3)),
    R = tibble(Q = ..1, D = ..2)
  )})) %>%
  select(-Interactors)

# save data
write_rds(test_collections, paste0(here(), "/test_collections.rds"))
write_rds(interactors_all, paste0(here(), "/interactors_all.rds"))
write_rds(interactors_2pubs, paste0(here(), "/interactors_2pubs.rds"))

# # Make a test collection for each interaction with two or more publications.
# test_collections = bind_rows(
#   ints %>% filter(Num_pubs >= 2) %>% select(Query = ProteinA, Rel_doc = ProteinB),
#   ints %>% filter(Num_pubs >= 2) %>% select(Query = ProteinB, Rel_doc = ProteinA),
# ) %>%
#   inner_join(interactors, by = c("Query" = "Protein"))
# 
# # FUNCTION: make_test_collection(query, rel_doc, )
# 
# # Data to make test collections.
# # Query: query protein
# # Rel_documents: proteins that interact with with the query as reported by at least two publications.
# # Interactors: proteins that interact with the query according to some evidence.
# test_collections <- bind_rows(
#   ints %>% filter(Num_pubs >= 2) %>% select(Query = ProteinA, Rel_doc = ProteinB),
#   ints %>% filter(Num_pubs >= 2) %>% select(Query = ProteinB, ProteinB = ProteinA)
# ) %>%
#   nest(Rel_documents = -ProteinA) %>%
#   mutate(Rel_documents = map(Rel_documents, ~ .x$ProteinB %>% unique)) %>%
#   select(Query = ProteinA, Rel_documents) %>%
#   inner_join(interactors, by = c("Query" = "Protein"))
# 
# # make_test_collections
# # Arguments:
# #   1. query - the query of the test collection
# #   2. rel_docs - set from which relevant documents are sampled
# #   3. docs <- set form which non-relevant documents are sampled
# #   4. interactors <- all proteins that interact with the query
# #   5. num_rel: number of relevant documents
# #   6. num: number of documents
# make_test_collections <- function(query, rel_docs, docs, interactors, num_rel, num) {
#   
#   # number of test collections that will be made
#   num_tc <- floor(length(rel_docs) / num_rel)
#   
#   if (num_tc < 1) 
#     return(NA)
#   
#   # initialize list of test collections
#   tcs <- vector(mode = "list", length = num_tc)
#   
#   # permute the relevant documents so that we can randomly distribute them in the test collections
#   rel_docs <- sample(rel_docs, size = length(rel_docs)) 
#   
#   # Make the test collections
#   for (i in 1:num_tc) {
#     
#     # range of relevant documents
#     end <- num_rel * i
#     start <- end - num_rel + 1
#     
#     # relevant documents of test collection
#     D_rel <- rel_docs[start:end]
#     
#     # relevant documents that are not chosen for test collection
#     rel_nc <- rel_docs[-(start:end)]
#     
#     # Sample non-relevant documents for test collection
#     D_nrel <- vector(mode = "character", length = num - num_rel)
#     j <- 1
#     while ( j <= (num - num_rel) ) {
#       
#       d <- sample(docs, 1)
#       
#       if ( !(d %in% D_nrel) && !(d %in% interactors) ) {
#         D_nrel[j] <- d
#         j <- j + 1
#       }
#     }
#     
#     # Sample some proteins that do not interact with the query
#     prots <- sample(docs) # shuffle proteins
#     # initialize proteins that do not interact
#     prots_ni <- vector(mode = "character", length = length(rel_nc))
#     j <- 1 # counter for prots
#     k <- 1 # counter for prots_ni
#     while( k <= length(rel_nc) ) {
#       
#       # sample protein
#       prot <- prots[j]
#       
#       # if not interactor of query
#       if ( !(prot %in% interactors) ) {
#         prots_ni[k] <- prot
#         k <- k + 1
#       }
#       
#       j <- j + 1
#     }
#     
#     # Form test collection.
#     # Pos: relevant proteins that were not chosen.
#     # Neg: Random sample of proteins that are not interactors (same size as Pos)
#     tcs[[i]] <- list(
#       Q = query,
#       D = c(D_rel, D_nrel),
#       R = tibble(Query = query, Document = D_rel),
#       Pos = rel_nc,
#       Neg = prots_ni
#     )
#   }
#   
#   return(tcs)
# }

# # 36 relevant documents
# query <- "O14920"
# rel_docs <- (test_collections %>% filter(Query == query))$Rel_documents[[1]]
# docs <- mparams$Protein
# interactors <- (test_collections %>% filter(Query == query))$Interactors[[1]]
# num_rel <- 1
# num <- 20
# tcs <- make_test_collections(query, rel_docs, docs, interactors, num_rel, num)
# 
# # 1 relevant document
# query <- "Q9UBE0"
# rel_docs <- (test_collections %>% filter(Query == query))$Rel_documents[[1]]
# docs <- mparams$Protein
# interactors <- (test_collections %>% filter(Query == query))$Interactors[[1]]
# num_rel <- 1
# num <- 20
# tcs1 <- make_test_collections(query, rel_docs, docs, interactors, num_rel, num)
# 
# # 2 relevant documents
# query <- "P20333"
# rel_docs <- (test_collections %>% filter(Query == query))$Rel_documents[[1]]
# docs <- mparams$Protein
# interactors <- (test_collections %>% filter(Query == query))$Interactors[[1]]
# num_rel <- 1
# num <- 20
# tcs2 <- make_test_collections(query, rel_docs, docs, interactors, num_rel, num)
# 
# # 3 relevant documents
# query <- "O75496"
# rel_docs <- (test_collections %>% filter(Query == query))$Rel_documents[[1]]
# docs <- mparams$Protein
# interactors <- (test_collections %>% filter(Query == query))$Interactors[[1]]
# num_rel <- 1
# num <- 20
# tcs3 <- make_test_collections(query, rel_docs, docs, interactors, num_rel, num)
# 
# # 4 relevant documents
# query <- "Q13541"
# rel_docs <- (test_collections %>% filter(Query == query))$Rel_documents[[1]]
# docs <- mparams$Protein
# interactors <- (test_collections %>% filter(Query == query))$Interactors[[1]]
# num_rel <- 1
# num <- 20
# tcs4 <- make_test_collections(query, rel_docs, docs, interactors, num_rel, num)
# 
# # 5 relevant documents
# query <- "Q14653"
# rel_docs <- (test_collections %>% filter(Query == query))$Rel_documents[[1]]
# docs <- mparams$Protein
# interactors <- (test_collections %>% filter(Query == query))$Interactors[[1]]
# num_rel <- 1
# num <- 20
# tcs5 <- make_test_collections(query, rel_docs, docs, interactors, num_rel, num)

# # Test collections
# all_docs <- mparams$Protein
# set.seed(27)
# pb <- progress_bar$new(total = nrow(test_collections))
# test_collections <- test_collections %>%
#   mutate(TCs = pmap(list(Query, Rel_documents, Interactors), ~ {
#     pb$tick();
#     make_test_collections(
#       query = ..1,
#       rel_docs = ..2,
#       docs = all_docs,
#       interactors = ..3,
#       num_rel = 1,
#       num = 20
#     )
#   }))
# 
# write_rds(test_collections, paste0(here(), "/test_collections.rds"))

