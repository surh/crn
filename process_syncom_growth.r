# setwd("/home/sur/lab/exp/2025/today3")
library(tidyverse)

source("/home/sur/lab/src/crn/functions.r")

Dat <- bind_rows(read_single_experiment(od600_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/od600.tsv",
                                        timepoints_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/timepoints.tsv",
                                        mpn_file = NULL,
                                        syncoms_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/syncoms.tsv",
                                        type = "syncom",
                                        batch_name = "NS1"),
                 read_single_experiment(od600_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS2/od600.tsv",
                                        timepoints_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS2/timepoints.tsv",
                                        mpn_file = NULL,
                                        syncoms_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS2/syncoms.tsv",
                                        type = "syncom",
                                        batch_name = "NS2"),
                 read_single_experiment(od600_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS3/od600.tsv",
                                        timepoints_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS3/timepoints.tsv",
                                        mpn_file = NULL,
                                        syncoms_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS3/syncoms.tsv",
                                        type = "syncom",
                                        batch_name = "NS3"),
                 read_single_experiment(od600_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS4/od600.tsv",
                                        timepoints_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS4/timepoints.tsv",
                                        mpn_file = NULL,
                                        syncoms_file = "/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS4/syncoms.tsv",
                                        type = "syncom",
                                        batch_name = "NS4")
)



write_tsv(Dat,"pilot_syncom_growth_curves.tsv")



