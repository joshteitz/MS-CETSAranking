library(conflicted)
library(tidyverse)
library(readxl)
library(here)
library(drc)
library(progress)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

### Introduction ###################################################################################################
# This file contains functions for loading data from the "Supplementary Materials" of the paper
# "Thermal proximity coaggregation for system-wide profiling of protein complex dynamics in cells" by Tan et al..
# This paper is accessible from the following url: https://science.sciencemag.org/content/359/6380/1170.
#
# In order to use this file:
#   1. Make sure there is a folder called "data" in the same location as "load_data.R"
#   2. Download "Tables S1 to S27" from https://science.sciencemag.org/content/suppl/2018/02/07/science.aan0346.DC1
#   2. Unzip folder and place it into "data" folder.
#   3. Ensure the folder is named "TablesS1_to_S27" and contains two Excel files (.xlsx) named
#      "Tables_S1_to_S18.xlsx" and "Tables_S19_to_S27.xlsx".
#   ** If either Excel file is open, the code will not run.
####################################################################################################################

# Location of the two Excel spreadsheets.
LOC <- here("data", "TablesS1_to_S27")

# Load melting data from a cell line.
load_mdata <- function(cell_line = "K562") {
  file_sheet = switch(
    cell_line,
    "K562 lysate" = c("Tables_S1_to_S18.xlsx", "Table S6"),
    "K562" = c("Tables_S1_to_S18.xlsx", "Table S7"),
    "A375" = c("Tables_S19_to_S27.xlsx", "Table S19"),
    "HCT116" = c("Tables_S19_to_S27.xlsx", "Table S20"),
    "HEK293T" = c("Tables_S19_to_S27.xlsx", "Table S21"),
    "HL60" = c("Tables_S19_to_S27.xlsx", "Table S22"),
    "MCF7" = c("Tables_S19_to_S27.xlsx", "Table S23"),
    "mouse liver" = c("Tables_S19_to_S27.xlsx", "Table S24")
  )
  if (is.null(file_sheet)) {
    stop("Unrecognizable cell line! Please enter one of: \"K562 lysate\", \"K562\", \"A375\", \"HCT116\",
         \"HEK293T\", \"HL60\", \"MCF7\", or \"mouse liver\".")
  }
  mdata <- read_excel(paste0(LOC, "/", file_sheet[1]), sheet = file_sheet[2], skip = 2) %>%
    select(Protein = Accession, matches("^T\\d+$"))
  return(mdata)
}

# Load protein-protein interactions 
load_interactions <- function(cell_line = c("S2", "293T", "HCT116")) {
  
  if (cell_line[1] == "S2") {
    
    # file/sheet that contains interactions
    file_name <- "Tables_S1_to_S18.xlsx"
    sheet_name <- "Table S2"
    
    # load data
    ints <- read_excel(paste0(LOC, "/", file_name), sheet = sheet_name, skip = 2) %>%
      select(ProteinA = `Protein A`, ProteinB = `Protein B`, Num_pubs = `No. of Publications`)
    
    # make sure ProteinA < ProteinB and make sure there are no duplicate interactions.
    ints <- ints %>% mutate(Proteins = map2(ProteinA, ProteinB, ~ sort(c(.x, .y)))) %>%
      mutate(ProteinA = map_chr(Proteins, ~ .x[1])) %>%
      mutate(ProteinB = map_chr(Proteins, ~ .x[2])) %>%
      select(-Proteins) %>%
      distinct(ProteinA, ProteinB, .keep_all = T)
    
  } else if (cell_line[1] %in% c("293T", "HCT116")) {
    
    file_name = switch(
      cell_line[1],
      "293T" = "interactions_293T.rds",
      "HCT116" = "interactions_hct116.rds"
    )
    
    ints <- read_rds(here("data", "BioPlex", file_name))
    
  } else {
    stop("Invalid cell line entered")
  }
  
  return(ints)
}

# Load CORUM protein complex data
load_cdata <- function() {
  cdata <- read_excel(paste0(LOC, "/Tables_S1_to_S18.xlsx"), sheet = "Table S8", skip = 2) %>%
    select(Complex = Complex.id, Proteins = subunits..UniProt.IDs.) %>%
    mutate(Complex = Complex %>% as.integer) %>%
    mutate(Proteins = strsplit(Proteins, ","))
  cdata <- map2_dfr(cdata$Complex, cdata$Proteins, ~ tibble(Complex = .x, Protein = .y))
  return(cdata)
}

# For each protein's melting data, fit a three parameter log-logistic function with upper limit equal to 1.
# Then extract the function's three parameters.
param_mdata <- function(mdata) {
  pb <- progress_bar$new(total = nrow(mdata))
  mparams <- convert_long(mdata) %>%
    group_by(Protein) %>%
    nest() %>%
    mutate(Mod = map(data, ~ {pb$tick(); get_drm(.x)})) %>%
    select(-data) %>%
    ungroup() %>%
    mutate(Param_b = map_dbl(Mod, ~ .x$coefficients[1])) %>%
    mutate(Param_c = map_dbl(Mod, ~ .x$coefficients[2])) %>%
    mutate(Param_e = map_dbl(Mod, ~ .x$coefficients[3])) %>%
    select(-Mod)
  return(mparams)
}

# Make specified columns unit variance.
make_unit_var <- function(mdata, numeric_cols = NULL) {
  scale1 <- function(x) {
    x_sd <- sd(x)
    x_mean <- mean(x)
    return((x - x_mean) / x_sd)
  }
  mdata_unit_var <- 
    mdata %>%
    mutate(across(any_of({{numeric_cols}}), scale1))
  return(mdata_unit_var)
}

# Convert melting data to long format
convert_long <- function(mdata) {
  mdata_long <- mdata %>%
    pivot_longer(
      cols = matches("T\\d+$"), 
      names_to = "Temp", 
      values_to = "Sol", 
      names_prefix = "T", 
      names_transform = list(Temp = as.integer)
    )
  return(mdata_long)
}

# Wrapper around `drm` function that fails gracefully when a dose response model cannot be made.
get_drm <- function(df) {
  tryCatch(
    drm(Sol ~ Temp, data = df, fct = LL.3u()),
    error = function(e) {
      list(coefficients = c(NA, NA, NA))
    }
  )
}

# FUNCTION: load proteins that come from a particular subcellular location
load_prots_loc <- function(loc) {
  
  file_name = switch(
    loc,
    "nucleoplasm" = "subcell_loc_nucleoplasm.tsv",
    "nuclear-membrane" = "subcell_loc_nuclear_membrane.tsv",
    "nucleoli" = "subcell_loc_nucleoli.tsv",
    "actin-filaments" = "subcell_loc_actin_filaments.tsv",
    "intermediate-filaments" = "subcell_loc_intermediate_filaments.tsv",
    "centrosome" = "subcell_loc_centrosome.tsv",
    "microtubules" = "subcell_loc_microtubules.tsv",
    "cytosol" = "subcell_loc_cytosol.tsv",
    "mitochondria" = "subcell_loc_mitochondria.tsv",
    "endoplasmic-reticulum" = "subcell_loc_endoplasmic_reticulum.tsv",
    "vesicles" = "subcell_loc_vesicles.tsv",
    "golgi-apparatus" = "subcell_loc_golgi_apparatus.tsv",
    "plasma-membrane" = "subcell_loc_plasma_membrane.tsv"
  )
  
  if (is.null(file_name)) {
    stop("Invalid subcellular location given to `loc` argument!")
  }
  
  prots <- read_tsv(here("data", "HumanProteinAtlas", file_name)) %>% pull(Uniprot)
  # prots <- read_tsv(here("data", "HumanProteinAtlas", file_name))
  return(prots)
}

# # FUNCTION: load bait-prey data
# #   fn: file name
# load_bait_prey <- function(fn) {
#   
#   file_name = switch(
#     fn,
#     "Ewings07" = "Ewings07.xls",
#     "BioPlex293T" = "BioPlexBaitPrey293T.tsv"
#   )
#   
#   if (is.null(file_name)) {
#     stop("Invalid file name given!")
#   }
#   
#   if (fn == "Ewings07") {
#     bp <- read_excel(here("data", file_name))
#   } else if (fn == "BioPlex293T") {
#     bp <- read_tsv(here("data", file_name))  
#   }else {
#     bp <- read_excel(here("data", file_name))
#   }
#   
#   return(bp)
# }

# FUNCTION: load BioPlex 3.0 data
#   fn: file name
load_bioplex <- function(fn) {
  
  file_name = switch(
    fn,
    "bait_prey_293T" = "bait_prey_293T.tsv",
    "bait_prey_hct116" = "bait_prey_hct116.tsv",
    "interactions_293T" = "interactions_293T.tsv",
    "interactions_hct116" = "interactions_hct116.tsv"
  )
  
  if (is.null(file_name)) {
    stop("Invalid file name given!")
  }
  
  if (fn %in% c("interactions_293T", "interactions_hct116")) {
    d <- read_tsv(here("data", "BioPlex", file_name)) %>%
      select(GeneA, SymbolA, GeneB, SymbolB, pNI, pInt)
  } else if (fn %in% c("bait_prey_293T", "bait_prey_hct116")) {
    d <- read_tsv(here("data", "BioPlex", file_name)) %>%
      select(bait_symbol, bait_geneid, symbol, gene_id)
  }
  
  return(d)
}

# FUNCTION: load bait/prey data
load_bait_prey <- function(cell_line) {
  
  file_name = switch(
    cell_line,
    "293T" = "bait_prey_293T.rds",
    "HCT116" = "bait_prey_hct116.rds"
  )
  
  if (is.null(file_name)) stop("Invalid file name given!")
  
  bp <- read_rds(here("data", "BioPlex", file_name))
  return(bp)
}
