# setwd("/home/sur/lab/exp/2025/today3")
library(tidyverse)

source("/home/sur/lab/src/crn/functions.r")


Dat <- bind_rows(read_single_experiment(od600_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/SHP1/od600.tsv",
                       timepoints_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/SHP1/timepoints.tsv",
                       mpn_file = NULL,
                       syncoms_file = NULL,
                       type = "strain",
                       batch_name = "SHP1"),
                read_single_experiment(od600_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/NS1/od600.tsv",
                       timepoints_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/NS1/timepoints.tsv",
                       mpn_file = NULL,
                       syncoms_file = NULL,
                       type = "strain",
                       batch_name = "NS1"),
                read_single_experiment(od600_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/ML1/od600.tsv",
                       timepoints_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/ML1/timepoints.tsv",
                       mpn_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/ML1/mpn.tsv",
                       syncoms_file = NULL,
                       type = "strain",
                       batch_name = "ML1")
)

strains_in_syncom <- c(
  "ST00046",
  "ST00154",
  "ST00101",
  "ST00109",
  "ST00042",
  "ST00060",
  "ST00094",
  "ST00110",
  "ST00164",
  "ST00143"
)

#' We keep only the strains in the syncom, and remove the ones that have
#' problems with the batch, using a majority rule
Dat.filtered <- Dat %>%
  filter(strain %in% strains_in_syncom) %>%
  filter(!(strain == "ST00060" & batch == "ML1")) %>%
  filter(!(strain == "ST00109" & batch == "SHP1")) %>%
  filter(!(strain == "ST00110" & batch == "SHP1")) 


# Count number of replicates per strain
Dat.filtered %>%
  filter(timepoint == "t_0") %>%
  group_by(strain, temp) %>%
  summarise(n_replicates = n()) %>%
  arrange(strain, temp) 


write_tsv(Dat,"pilot_strain_growth_curves.tsv")
write_tsv(Dat.filtered,"pilot_strain_growth_curves_filtered.tsv")



