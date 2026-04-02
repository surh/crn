setwd("/home/sur/lab/exp/2026/today/")
library(tidyverse)
library(brms)
box::use(./functions/crn)


#' # Simulate some simple data
#' A bit cumbersome, but the line below simulate data for 4 **id**'s
#' (could be strains or syncoms) in 3 temps, 2 reps, and 4 batches (1 rep
#' of 2 id's per batch). A switch in id's per batch is simulated.
#' Overall variability is modelled (includes measurement error) and
#' batch-to-batch variations is simulated as additive 4-times smaller
#' than replicate variability (in reality it could be opposite, more batch
#' than biological variability)
set.seed(1234567)
var_rep <- 0.2
var_batch <- 0.05
b_batch <- rnorm(n = 4, mean = 0, sd = sqrt(var_batch))
names(b_batch) <- c("A","B","C","D")
Dat <- bind_rows(
  # First rep
  tibble(id = "ST01",
         temp = c(30,37,42),
         AUC = c(4,4,4),
         batch = "A"),
  tibble(id = "ST02",
         temp = c(30,37,42),
         AUC = c(4,4.5,5),
         batch = "A"),
  tibble(id = "ST03",
         temp = c(30,37,42),
         AUC = c(4,4.5,3.7),
         batch = "B"),
  
  tibble(id = "ST04",
         temp = c(30,37,42),
         AUC = c(4,4.1,3.2),
         batch = "B"),
  
  # Second rep with flipped batches
  tibble(id = "ST01",
         temp = c(30,37,42),
         AUC = c(4,4,4),
         batch = "C"),
  tibble(id = "ST02",
         temp = c(30,37,42),
         AUC = c(4,4.5,5),
         batch = "D"),
  tibble(id = "ST03",
         temp = c(30,37,42),
         AUC = c(4,4.5,3.7),
         batch = "C"),
  tibble(id = "ST04",
         temp = c(30,37,42),
         AUC = c(4,4.1,3.2),
         batch = "D")
) %>%
  mutate(AUC = AUC + rnorm(n = length(id), mean = 0, sd = var_rep)) %>%
  mutate(AUC = AUC + b_batch[batch])

#' Contigency table of samples
ftable(id + batch ~ temp, data = Dat)

#' Plot of the resulting *reaction norms*.
p1 <- Dat %>%
  ggplot(aes(x = temp, y = AUC,
             group = interaction(id, batch, sep = "_", drop = TRUE))) +
  geom_line(aes(col = id), linewidth = 3) +
  theme_classic()
p1

#' # Reaction norm
#' First it is convenient to convert temp to -1, 0, 1 for low, med and high

Dat <- Dat %>%
  mutate(temp = replace(temp, temp == 30, -1)) %>%
  mutate(temp = replace(temp, temp == 37, 0)) %>%
  mutate(temp = replace(temp, temp == 42, 1))
Dat


mp0 <- brm(AUC ~ 1 + temp + temp^2 + (1 | id) + (1 | batch),
           data = Dat,
           chains = 4, iter = 4000, warmup = 3000, threads = 4,
           control = list(adapt_delta = 0.95))
summary(mp0)



mp1 <- brm(AUC ~ 1 + (1 | temp + temp^2) + (1 | id) + (1 | batch),
           data = Dat,
           chains = 4, iter = 4000, warmup = 3000)

summary(mp1)


mp2 <- brm(AUC ~ 1 + temp + (1 + temp | id), data = Dat)
m1.brms

















