library(tidyverse)
library(lme4)
setwd("asanchez_lab")
# source("../functions.r")

args <- list()
args$outdir <- "/home/sur/micropopgen/exp/2024/today/"

Dat <- read_tsv("data/pseudo.txt")
Dat


# First pass
Dat <- Dat %>%
    filter(wavelength == 600) %>%
    mutate(dilution_factor_sq = dilution_factor^2) %>%
    mutate(dilution_factor_f = factor(dilution_factor))



p1 <- ggplot(Dat, aes(
    x = 1/dilution_factor, y = absorbance,
    group = interaction(community))
) +
    # geom_line(alpha = 0.3) +
    # geom_point() +
    geom_smooth(method = "lm", se = FALSE, alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    labs(
        title = "Absorbance vs Dilution Factor",
        x = "Dilution Factor",
        y = "Absorbance"
    )

p1
ggsave(file.path(args$outdir, "absorbance_vs_dilution_factor_crn.png"), p1, width = 6, height = 4)




Crn.full <- file.path(args$outdir, "full_data_crns.tsv") %>%
    read_tsv()

Crn.rand <- file.path(args$outdir, "random_draws_crns.tsv") %>%
    read_tsv()

Crn.f2 <- file.path(args$outdir, "f2like_crns.tsv") %>%
    read_tsv(col_types = cols(
        expid = col_character()
    ))


p1 <- Crn.rand %>%
    mutate(strategy = "random") %>%
    bind_rows(
        Crn.f2 %>%
            mutate(strategy = "f2like")
    ) %>%
    ggplot(aes(x = V_Tot)) +
    facet_wrap(~strategy, scales = "free_y") +
    geom_histogram(bins = 20) +
    # geom_density(n = 256) +
    geom_vline(xintercept = Crn.full$V_Tot, col = "red", size = 2) +
    theme_classic()
p1
ggsave(file.path(args$outdir, "V_Tot_crn.png"), p1, width = 6, height = 4)



p1 <- Crn.rand %>%
    mutate(strategy = "random") %>%
    bind_rows(
        Crn.f2 %>%
            mutate(strategy = "f2like")
    ) %>%
    ggplot(aes(x = V_Phen)) +
    facet_wrap(~strategy, scales = "free_y") +
    geom_histogram(bins = 20) +
    # geom_density(n = 256) +
    geom_vline(xintercept = Crn.full$V_Phen, col = "red", size = 2) +
    theme_classic()
p1
ggsave(file.path(args$outdir, "V_Phen_crn.png"), p1, width = 6, height = 4)



p1 <- Crn.rand %>%
    mutate(strategy = "random") %>%
    bind_rows(
        Crn.f2 %>%
            mutate(strategy = "f2like")
    ) %>%
    ggplot(aes(x = V_Gen)) +
    facet_wrap(~strategy, scales = "free_y") +
    geom_histogram(bins = 10) +
    # geom_density(n = 256) +
    geom_vline(xintercept = Crn.full$V_Gen, col = "red", size = 2) +
    theme_classic()
p1
ggsave(file.path(args$outdir, "V_Gen_crn.png"), p1, width = 6, height = 4)
    
    

p1 <- Crn.rand %>%
    mutate(strategy = "random") %>%
    bind_rows(
        Crn.f2 %>%
            mutate(strategy = "f2like")
    ) %>%
    ggplot(aes(x = V_Res)) +
    facet_wrap(~strategy, scales = "free_y") +
    geom_histogram(bins = 10) +
    # geom_density(n = 256) +
    geom_vline(xintercept = Crn.full$V_Res, col = "red", size = 2) +
    theme_classic()
p1
ggsave(file.path(args$outdir, "V_Res_crn.png"), p1, width = 6, height = 4)


p1 <- Crn.rand %>%
    mutate(strategy = "random") %>%
    bind_rows(
        Crn.f2 %>%
            mutate(strategy = "f2like")
    ) %>%
    ggplot(aes(x = V_Plas)) +
    facet_wrap(~strategy, scales = "free_y") +
    geom_histogram(bins = 10) +
    # geom_density(n = 256) +
    geom_vline(xintercept = Crn.full$V_Plas, col = "red", size = 2) +
    theme_classic()
p1
ggsave(file.path(args$outdir, "V_Plas_crn.png"), p1, width = 6, height = 4)


