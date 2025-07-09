library(tidyverse)
library(coda)

box::use(./fit_growth/functions)


args <- list()
args$indir <- "/Users/sur/lab/exp/2025/today3/single_strain_bytemp_logistic/"
infiles <- list.files(args$indir, pattern = "_logistic_fit.tsv", full.names = TRUE)

Res <- tibble()
for(f in infiles){
    # f <- infiles[1]
    strain <- basename(f) %>% str_remove("_logistic_fit.tsv") %>% str_split("_")
    temp <- strain[[1]][2]
    strain <- strain[[1]][1]
    cat(strain, "\t", temp, "\n")

    Dat <- read_tsv(f)

    res <- bind_rows(functions$calculate_stats(Dat = Dat, column = "r"),
        functions$calculate_stats(Dat = Dat, column = "K")) %>%
        mutate(strain = strain) %>%
        mutate(temp = temp) %>%
        select(strain, temp, everything())

    Res <- bind_rows(Res, res)

}
Res

Dat
Dat %>%
    select(K, chain, iteration) %>%
    pivot_wider(names_from = chain, values_from = K) %>%
    select(-iteration) %>%
    as.matrix() %>%
    rstan::Rhat()



p1 <- Res %>%
    filter(parameter == "r") %>%
    ggplot(aes(x = temp, y = mean)) +
    facet_grid(. ~ strain, scales = "free_y") +
    geom_linerange(aes(ymin = hpd05, ymax = hpd95), linewidth = 1, color = "darkgrey") +
    geom_linerange(aes(ymin = hpd10, ymax = hpd90), linewidth = 3, color = "darkgrey") +
    geom_point(size = 6, color = "#dd1c77") +
    geom_point(aes(y = median), size = 5, color = "#c994c7") +
    labs(y = "Growth rate (1/h)", x = "Temperature (ºC)") +
    theme_classic()
p1


p1 <- Res %>%
    filter(parameter == "K") %>%
    ggplot(aes(x = mean, y = temp)) +
    facet_wrap(~strain, scales = "fixed") +
    geom_linerange(aes(xmin = hpd05, xmax = hpd95), size = 1, color = "darkgrey") +
    geom_linerange(aes(xmin = hpd10, xmax = hpd90), size = 3, color = "darkgrey") +
    geom_point(size = 6, color = "#dd1c77") +
    geom_point(aes(x = median), size = 5, color = "#c994c7") +
    labs(x = "Carrying capacity (OD600)", y = "Strain") +
    theme_classic()
p1
