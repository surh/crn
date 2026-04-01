setwd("/home/sur/lab/exp/2026/today/")
library(tidyverse)
box::use(../../../src/crn/functions/crn)


#' # Simulate some simple data

expand_grid(strain = letters[1:10],
            temp = c(-1,0,1),
            rep = 1:2) %>%
  mutate(batch = sample(LETTERS[1:4], size = length(strain), replace = TRUE)) %>%
  mutate(AUC = rnorm(2) + 1.5 * temp + 3 * temp^2) %>%
  ggplot(aes(x = temp, y = AUC,
             colour = interaction(strain, rep, batch, sep = "_"))) +
  geom_line()




