library(tidyverse)
library(coda)

calculate_stats <- function(Dat, column = "r"){

    par.mcmc <- as.mcmc(Dat[[column]])
    hpd <- HPDinterval(par.mcmc, prob = 0.9)
    hpd05 <- hpd[1]
    hpd95 <- hpd[2]
    hpd <- HPDinterval(par.mcmc, prob = 0.8)
    hpd10 <- hpd[1]
    hpd90 <- hpd[2]

    tibble(parameter = column,
           mean = mean(as.numeric(Dat[[column]]), na.rm = TRUE),
           sd = sd(as.numeric(Dat[[column]]), na.rm = TRUE),
           median = median(as.numeric(Dat[[column]]), na.rm = TRUE),
           q05 = quantile(as.numeric(Dat[[column]]), 0.05, na.rm = TRUE),
           q10 = quantile(as.numeric(Dat[[column]]), 0.1, na.rm = TRUE),
           q90 = quantile(as.numeric(Dat[[column]]), 0.9, na.rm = TRUE),
           q95 = quantile(as.numeric(Dat[[column]]), 0.95, na.rm = TRUE),
           hpd05 = hpd05,
           hpd10 = hpd10,
           hpd90 = hpd90,
           hpd95 = hpd95
    )
}

args <- list()
indir <- "/Users/sur/lab/exp/2025/today3/single_strain_logistic/"
infiles <- list.files(indir, pattern = "_logistic_fit.tsv", full.names = TRUE)

Res <- tibble()
for(f in infiles){
    # f <- infiles[1]
    strain <- basename(f) %>% str_remove("_logistic_fit.tsv")
    cat(strain, "\n")

    Dat <- read_tsv(f)

    res <- bind_rows(calculate_stats(Dat = Dat, column = "r"),
        calculate_stats(Dat = Dat, column = "K")) %>%
        mutate(strain = strain) %>%
        select(strain, everything())

    Res <- bind_rows(Res, res)

}
Res


p1 <- Res %>%
    filter(parameter == "r") %>%
    ggplot(aes(x = mean, y = strain)) +
    geom_errorbarh(aes(xmin = hpd05, xmax = hpd95), size = 1, color = "darkgrey") +
    geom_linerange(aes(xmin = hpd10, xmax = hpd90), size = 3, color = "darkgrey") +
    geom_point(size = 6, color = "#dd1c77") +
    geom_point(aes(x = median), size = 5, color = "#c994c7") +
    labs(x = "Growth rate (1/h)", y = "Strain") +
    theme_classic()
p1

p1 <- Res %>%
    filter(parameter == "K") %>%
    ggplot(aes(x = mean, y = strain)) +
    geom_errorbarh(aes(xmin = hpd05, xmax = hpd95), size = 1, color = "darkgrey") +
    geom_linerange(aes(xmin = hpd10, xmax = hpd90), size = 3, color = "darkgrey") +
    geom_point(size = 6, color = "#dd1c77") +
    geom_point(aes(x = median), size = 5, color = "#c994c7") +
    labs(x = "Carrying capacity (OD600)", y = "Strain") +
    theme_classic()
p1
