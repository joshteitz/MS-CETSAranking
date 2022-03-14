library(here)

source(here("load_data.R"))

# FUNCTION: make test collection. R is created by filtering all_rel_judgs
make_tc <- function(query, documents, all_rel_judgs) {
  
  # subcellular location(s) of Q, possibly NULL
  loc <- prots_loc %>% filter(Protein == query) %>% pull(Loc) %>% unlist()
  
  # D consists only of documents that do not belong to a different subcellular location than
  # the query.
  docs <- tibble(Protein = documents) %>% 
    left_join(prots_loc, by = "Protein") %>%
    mutate(Loc_compatible = map_lgl(Loc, ~ {
      if (is.null(loc) || is.null(.x))
        TRUE
      else if (intersect(loc, .x) %>% length > 0)
        TRUE
      else 
        FALSE
    })) %>%
    filter(Loc_compatible) %>%
    pull(Protein)
  
  R <- all_rel_judgs %>% filter((Q == query) & (D %in% docs))
  
  if (nrow(R) > 0) {
    tc <- list(Q = query, D = docs, R = R)
  } else {
    tc <- NULL
  }
  
  return(tc)
}

# FUNCTION: make all test collections
#   cell_bp: cell line of bait/prey data
#   cell_ints: cell line of interactions
make_all_tcs <- function(
  cell_bp = c("293T", "HCT116"),
  cell_ints = c("S2", "293T", "HCT116")
) {
  
  # # remove proteins without parametric melting data
  # prots_loc <- prots_loc %>% semi_join(mparams, by = "Protein")
  
  # load bait/prey data
  bp <- load_bait_prey(cell_bp[1])
  
  # # remove bait/prey data without parametric melting data
  # bp <- bp %>% 
  #   semi_join(mparams, by = c("Bait" = "Protein")) %>%
  #   semi_join(mparams, by = c("Prey" = "Protein"))
  
  # remove rows where bait and prey are the same
  bp <- bp %>% filter(Bait != Prey) 
  
  # make sure all rows are distinct
  bp <- bp %>% distinct(Bait, Prey)
  
  # interactions
  if (cell_ints[1] == "S2") {
    ints <- load_interactions(cell_ints[1]) %>% filter(Num_pubs >= 2)
  } else {
    ints <- load_interactions(cell_ints[1]) %>% filter(pInt > .95)
  }
  
  # # Remove each interaction where at least one protein does not have parametric melting data
  # ints <- ints %>%
  #   semi_join(mparams, by = c("ProteinA" = "Protein")) %>%
  #   semi_join(mparams, by = c("ProteinB" = "Protein"))
  
  # relevance judgements from which R is drawn
  all_rel_judgs <- all_rel_judgs <- bind_rows(
    tibble(Q = ints$ProteinA, D = ints$ProteinB),
    tibble(Q = ints$ProteinB, D = ints$ProteinA)
  ) %>%
    distinct(Q, D)
  
  # make test collections
  print("Making a test collection for each bait...")
  pb <- progress_bar$new(total = bp$Bait %>% unique %>% length)
  tcs <- bp %>% nest(Preys = Prey) %>%
    mutate(Preys = map(Preys, ~ .x$Prey)) %>%
    mutate(TC = map2(Bait, Preys, ~ {pb$tick(); make_tc(.x, .y, all_rel_judgs)})) %>%
    drop_na(TC) %>%
    select(-Preys)
  
  return(tcs)
}

# Load protein location data
prots_loc <- bind_rows(
  tibble(Protein = load_prots_loc("nucleoplasm"), Location = "nucleoplasm"),
  tibble(Protein = load_prots_loc("nuclear-membrane"), Location = "nuclear-membrane"),
  tibble(Protein = load_prots_loc("nucleoli"), Location = "nucleoli"),
  tibble(Protein = load_prots_loc("actin-filaments"), Location = "actin-filaments"),
  tibble(
    Protein = load_prots_loc("intermediate-filaments"),
    Location = "intermediate-filaments"
  ),
  tibble(Protein = load_prots_loc("centrosome"), Location = "centrosome"),
  tibble(Protein = load_prots_loc("microtubules"), Location = "microtubules"),
  tibble(Protein = load_prots_loc("cytosol"), Location = "cytosol"),
  tibble(Protein = load_prots_loc("mitochondria"), Location = "mitochondria"),
  tibble(
    Protein = load_prots_loc("endoplasmic-reticulum"),
    Location = "endoplasmic-reticulum"
  ),
  tibble(Protein = load_prots_loc("vesicles"), Location = "vesicles"),
  tibble(Protein = load_prots_loc("golgi-apparatus"), Location = "golgi-apparatus"),
  tibble(Protein = load_prots_loc("plasma-membrane"), Location = "plasma-membrane"),
) %>%
  drop_na(Protein) %>%
  nest(Loc = Location) %>%
  mutate(Loc = map(Loc, ~ .x$Location %>% unique %>% sort))

# Bait/prey data: 293T
# Interactions: S2 with more than two citations
tcs1 <- make_all_tcs(cell_bp = "293T", cell_ints = "S2")
write_rds(tcs1, here("tcs", "tcs1.rds"))

# Bait/prey data: 293T
# Interactions: 293T
tcs2 <- make_all_tcs(cell_bp = "293T", cell_ints = "293T")
write_rds(tcs2, here("tcs", "tcs2.rds"))

# Bait/prey data: HCT116
# Interactions: S2 with more than two citations
tcs3 <- make_all_tcs(cell_bp = "HCT116", cell_ints = "S2")
write_rds(tcs3, here("tcs", "tcs3.rds"))

# Bait/prey data: HCT116
# Interactions: HCT116
tcs4 <- make_all_tcs(cell_bp = "HCT116", cell_ints = "HCT116")
write_rds(tcs4, here("tcs", "tcs4.rds"))
