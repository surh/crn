# setwd("/Users/sur/lab/exp/2024/today2")
library(tidyverse)
# library(lme4)
# library(brms)

#' Takes matrices produced by Natalia said and prepares tibbles for model inputs.
#' 
#' First we setup some parameters:
args <- list()
args$outdir <- "data_from_gore_for_models"
args$pheno_matrix <- "Matrix_OD.csv"
args$count_matrix <- "Matrix_abundance.csv"
args$metadata_table <- "Matrix_reference.csv"

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

#' Read data. Clean ID names to remove spaces
Pheno <- read_csv(args$pheno_matrix)
Pheno

Counts <- read_csv(args$count_matrix)
names(Counts) <- c("Species", names(Counts)[-1] %>% str_remove(" "))
Counts

#' Read metadata table
Meta <- read_csv(args$metadata_table)
Meta <- Meta %>%
  mutate(ID = str_remove(ID, " "))
Meta

#' Match phenotype data with metadata (i.e. add IDs)
Pheno <- Pheno %>%
  full_join(Meta, by = c("Day", "Rep", "Community", "S0", "Nutrient"))
Pheno

#' Check that no ID appears more than once
if (!all(table(Pheno$ID[!is.na(Pheno$ID)]) == 1)){
  stop("ERROR. At least one ID is repeated in Pheno table")
}

#' Check that all IDs in Meta show up in Pheno
if(length(setdiff(Meta$ID, Pheno$ID[!is.na(Pheno$ID)])) != 0){
  stop("ERROR. At least one 1 from Meta missing in Pheno")
}

#' Generate M IDs for Pheno measurements without IDs
Pheno$ID[is.na(Pheno$ID)] <- paste0("M", 
                                    (max(Pheno$ID[!is.na(Pheno$ID)] %>%
                                           str_remove("^M") %>%
                                           as.numeric()) + 1):nrow(Pheno))
Pheno

#' Recode nutrient as -1 (LN), 0 (MN), and 1 (HN)
Pheno <- Pheno %>%
  rename(Nutrient_orig = Nutrient) %>%
  mutate(Nutrient = NA) %>%
  mutate(Nutrient = replace(Nutrient, Nutrient_orig == "LN", -1)) %>%
  mutate(Nutrient = replace(Nutrient, Nutrient_orig == "MN", 0)) %>%
  mutate(Nutrient = replace(Nutrient, Nutrient_orig == "HN", 1)) %>%
  select(-Nutrient_orig)
Pheno


#' Temporary fix. Create new community with combination of Community and S0
Pheno <- Pheno %>%
  mutate(Com = interaction(Community, S0, sep = "_") %>% as.character()) %>%
  select(-Community) %>%
  rename(Community = Com)
Pheno

#' Get data for day 1 
Dat <- Pheno %>%
  filter(Day == 1)
Dat
outfile <- file.path(args$outdir, "day1_data.tsv")
write_tsv(Dat, outfile)


#' **To do**: instead of modeling OD directly, we should model \delta OD,
#' meaning the change in OD either from day 0, or from the previous day.

