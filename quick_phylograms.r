# setwd("syneco/exp/2024/today4/")

library(tidyverse)
library(AMOR)


Tab <- read_tsv("data/feature-table-closed.tsv", skip = 1) 
colnames(Tab)[1] <- "strain"
colnames(Tab) <- colnames(Tab) %>% str_replace("^N5", "NS")
Tab


Meta <- read_tsv("data/collection_sheet_combined.tsv")
Meta <- Meta %>% filter(label %in% colnames(Tab)) %>%
  arrange(match(Meta$label, colnames(Tab)))
Meta

Strains <- read_tsv("data/strains.tsv")
Strains$color_code <- c(rev(RColorBrewer::brewer.pal(n = 5, name = "Blues")),
                        RColorBrewer::brewer.pal(n = 5, name = "Reds"))


Dat <- create_dataset(as.matrix(Tab[,-1]), Meta %>% dplyr::select(label, everything()))

Tab.mat <- matrix(as.matrix(Tab[,-1]), nrow = nrow(Tab), dimnames = list(Tab$strain, colnames(Tab)[-1]))
Meta.df <- as.data.frame(Meta)
row.names(Meta.df) <- Meta$label

p1 <- phylogram(Tab = Tab.mat, Map = Meta.df,
          facet = "community",
          variable.name = "strain")

dat <- p1$data

p1 <- dat %>%
  filter(temp %in% c(NA, 28)) %>%
  mutate(community = factor(community, levels = paste0("R", 1:12))) %>%
  mutate(strain = factor(strain, levels = Strains$strain)) %>%
  ggplot(aes(x = day, y = Abundance)) +
  facet_wrap(community ~ ., scales = "free_x", nrow = 1) +
  geom_bar(aes(fill = strain), stat = "identity", position = "fill") +
  scale_fill_manual(values = Strains$color_code) +
  theme_classic()
ggsave("syncoms_28_rep1.png", width = 8, height = 3)

p1 <- dat %>%
  filter(temp %in% c(NA, 32)) %>%
  mutate(community = factor(community, levels = paste0("R", 1:12))) %>%
  mutate(strain = factor(strain, levels = Strains$strain)) %>%
  ggplot(aes(x = day, y = Abundance)) +
  facet_wrap(community ~ ., scales = "free_x", nrow = 1) +
  geom_bar(aes(fill = strain), stat = "identity", position = "fill") +
  scale_fill_manual(values = Strains$color_code) +
  theme_classic()
ggsave("syncoms_32_rep1.png", width = 8, height = 3)
