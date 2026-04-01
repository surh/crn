setwd("/home/sur/lab/exp/2026/today/")
library(tidyverse)
box::use(./functions/crn)


#' # Simulate some simple data
options(device = function(...) png(type = "cairo", ...))
expand_grid(strain = letters[1:10],
            temp = c(-1,0,1),
            rep = 1:2) %>%
  mutate(batch = sample(LETTERS[1:4], size = length(strain), replace = TRUE)) %>%
  mutate(AUC = rnorm(2) + 1.5 * temp + 3 * temp^2) %>%
  ggplot(aes(x = temp, y = AUC,
             colour = interaction(strain, rep, batch, sep = "_"))) +
  geom_line() +
  theme(text = element_text(family = "sans"))

# op <- options(device = "cairo")
expand_grid(bind_rows(tibble(strain = "ST01",
                 temp = c(30,37,42),
                 AUC = c(4,4,4)),
          tibble(strain = "ST02",
                 temp = c(30,37,42),
                 AUC = c(4,4.5,5)),
          tibble(strain = "ST03",
                 temp = c(30,37,42),
                 AUC = c(4,4.5,3.7)),
          
          tibble(strain = "ST04",
                 temp = c(30,37,42),
                 AUC = c(4,4.1,3.2))
), rep = 1:2) %>%
  mutate(batch = sample(LETTERS[1:3], size = length(strain), replace = TRUE))


tibble(x = 1:10,
       y = 1:10) %>%
  ggplot(aes(x = x, y = y)) +
  geom_point()
