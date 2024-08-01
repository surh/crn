# setwd("/Users/sur/lab/exp/2024/today5")
library(tidyverse)

# OD <- read_tsv("od600_strain_test.tsv")
# Timepoints <- read_tsv("od600_strain_test_timepoints.tsv")




read_single_experiment <- function(od600_file, timepoints_file, syncoms_file = NULL,
                                   mpn_file = NULL, type = "syncom",
                                   batch_name =  basename(dirname(od600_file))){
  
  # od600_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/od600.tsv"
  # timepoints_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/timepoints.tsv"
  # syncoms_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/syncoms.tsv"
  # mpn_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/mpn.tsv"
  # type <- "syncom"
  # batch_name <- "NS1"

  if(type == "syncom"){
    id_cols <- c("Community", "temp")
  }else if(type == "strain"){
    id_cols <- c("Strain", "temp")
  }
  
  
  #' Read and proces OD
  OD <- read_tsv(od600_file)
  Timepoints <- read_tsv(timepoints_file)
  
  Dat <- OD %>%
    pivot_longer(-id_cols,
                 names_to = "timepoint",
                 values_to = "OD600") %>%
    left_join(Timepoints, by = c("timepoint", "temp")) %>%
    mutate(OD600 = (OD600 - blank_od) * od_factor) %>%
    mutate(OD600 = replace(OD600, timepoint == "t_0", 0.001)) 
  
  
  #' Add MPN if available
  if(!is.null(mpn_file) && file.exists(mpn_file)){
    MPN <- read_tsv(mpn_file, na = c("", "NA", "Error"))
    
    Dat <- MPN %>%
      pivot_longer(-id_cols,
                   names_to = "timepoint",
                   values_to = "MPN") %>%
      right_join(Dat, by = c(id_cols, "timepoint")) 
  }
  
  if(!is.null(batch_name)){
    Dat[["batch"]] <- batch_name
  }
  
  
  return(Dat)
  
}

#' Read all syncom data
Syncoms <- bind_rows(read_single_experiment(od600_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/od600.tsv",
                                            timepoints_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/timepoints.tsv",
                                            mpn_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/mpn.tsv",
                                            syncoms_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/syncoms.tsv",
                                            type = "syncom",
                                            batch_name = "NS1"),
                     read_single_experiment(od600_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS2/od600.tsv",
                                            timepoints_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS2/timepoints.tsv",
                                            mpn_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS2/mpn.tsv",
                                            syncoms_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS2/syncoms.tsv",
                                            type = "syncom",
                                            batch_name = "NS2"),
                     read_single_experiment(od600_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS3/od600.tsv",
                                            timepoints_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS3/timepoints.tsv",
                                            mpn_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS3/mpn.tsv",
                                            syncoms_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS3/syncoms.tsv",
                                            type = "syncom",
                                            batch_name = "NS3"),
                     read_single_experiment(od600_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS4/od600.tsv",
                                            timepoints_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS4/timepoints.tsv",
                                            mpn_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS4/mpn.tsv",
                                            syncoms_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS4/syncoms.tsv",
                                            type = "syncom",
                                            batch_name = "NS4"))

Syncoms
  







#' Remove t_2 which seems to have overestimated OD600
p1 <- Dat %>%
  # filter(timepoint != "t_2") %>% 
  ggplot(aes(x = total_time_h, y = OD600, group = interaction(Community, temp, sep = "_"))) +
  facet_wrap(~ Community) +
  geom_vline(xintercept = c(24, 48)) +
  geom_line(aes(col = as.character(temp))) +
  scale_color_manual(values = c("blue", "red")) +
  theme_classic()
p1
ggsave("od_plots.svg", width = 6, height = 4)


#' It seems like at dillution we are reaching OD
#' in the order of 1e3, and reaching same growth
#' than when starting with 1e2. Thus maybe it is worth starting
#' at 
OD %>%
  select(Community, temp, t_1) %>%
  mutate(t_1 = t_1 / 100) %>%
  print(n = 100)


#' Check changes in growth performance
dat <- Dat %>%
  group_by(Community, timepoint) %>%
  summarise(max = temp[which.max(OD600)],
            equal = OD600[temp == 28] == OD600[temp == 32],
            .groups = 'drop') %>%
  filter(!equal) %>%
  select(-equal) %>%
  filter(timepoint %in% c("t_1", "t_2")) %>%
  pivot_wider(names_from = "timepoint", values_from = "max")
dat

p1 <- dat %>%
  pivot_longer(-Community, names_to = "timepoint", values_to = "temp") %>%
  ggplot(aes(x = timepoint)) +
  geom_bar(aes(fill = as.character(temp))) +
  scale_fill_manual(values = c("blue", "red")) +
  ggtitle(label = "Temperature with highest OD600") +
  theme_classic()
p1
ggsave("winning_temp_bar.svg", width = 4, height = 4)

#' Interestingly the group high_slow is the largest group. In
#' All cases when there is a flip, there is a reduction in the OD
#' between t_3 and t_4 for the "fast" condtion. 
#' In fact the only case where both conditions increase their OD is ST00109
# dat %>%
#   mutate(class = ifelse(t_3 == 28 & t_4 == 28, "low_temp", NA)) %>%
#   mutate(class = replace(class, t_3 == 28 & t_4 == 32, "high_slow")) %>%
#   mutate(class = replace(class, t_3 == 32 & t_4 == 28, "low_slow")) %>%
#   mutate(class = replace(class, t_3 == 32 & t_4 == 32, "high_temp")) %>%
#   mutate(class = factor(class, levels = c("low_temp", "high_slow", "low_slow", "high_temp"))) %>%
#   arrange(class) %>% print(n = 100)
#' This suggests that increasing temperature changes the performance of strains
#' but only after passage?. Or that at 32ºC these strains can recover better from
#' death phase?



#' Interestingly for the time between t_4 t_5 there is a completel even
#' distribution of groups.
# dat %>%
#   mutate(class = ifelse(t_4 == 28 & t_5 == 28, "low_temp", NA)) %>%
#   mutate(class = replace(class, t_4 == 28 & t_5 == 32, "high_slow")) %>%
#   mutate(class = replace(class, t_4 == 32 & t_5 == 28, "low_slow")) %>%
#   mutate(class = replace(class, t_4 == 32 & t_5 == 32, "high_temp")) %>%
#   mutate(class = factor(class, levels = c("low_temp", "high_slow", "low_slow", "high_temp"))) %>%
#   arrange(class) %>% print(n = 100)

Dat

Dat.mpn  <- tibble(Community="R1",
                   temp = c(28,32, 28,32,28,32),
                   MPN = c(3465487, 3465487, 1.25e9, 2.31e9, 1.79e9, 2.037e9),
                   total_time_h = c(0,0,24,24,48,48))

p1 <- Dat.mpn %>%
  # filter(timepoint != "t_2") %>% 
  ggplot(aes(x = total_time_h, y = MPN, group = interaction(Community, temp, sep = "_"))) +
  facet_wrap(~ Community) +
  geom_vline(xintercept = c(24, 48)) +
  geom_line(aes(col = as.character(temp))) +
  scale_color_manual(values = c("blue", "red")) +
  # scale_y_log10() + 
  theme_classic()
p1

ggsave("od_plots.svg", width = 6, height = 4)
