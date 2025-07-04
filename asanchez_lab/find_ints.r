

library(tidyverse)
library(lme4)
setwd("asanchez_lab")
source("../functions.r")

args <- list()
args$outdir <- "/home/sur/micropopgen/exp/2024/today/"

Dat <- read_tsv("data/pseudo.txt")
Dat


# First pass
Dat <- Dat %>%
    filter(wavelength == 600) %>%
    mutate(dilution_factor_sq = dilution_factor^2) %>%
    mutate(dilution_factor_f = factor(dilution_factor))




# Now for the f2-like design

Crns <- combn(1:8, 4) %>%
    as_tibble() %>%
    as.list() %>%
    map_dfr(function(red_strains) {
        # red_strains <- c(1, 3, 5, 7)
        blue_strains <- setdiff(1:8, red_strains)

        # Make all possible orders of blue strains
        blue_permutations <- gtools::permutations(
            n = length(blue_strains),
            r = length(blue_strains), v = blue_strains
        )

        Res.per.red <- blue_permutations %>%
            t() %>%
            as_tibble() %>%
            as.list() %>%
            map_dfr(function(blue_perm) {
                # blue_perm <- c(4, 6, 2, 8)
                # blue_perm
                # red_strains


                red_com <- rep(0, 8)
                red_com[red_strains] <- 1
                # Select communities for the design
                Coms <- paste0(red_com, collapse = "")
                for (i in 1:4) {
                    coms <- apply(combn(red_strains, i), 2, function(ii, red_com, blue_perm, red_strains) {
                        new_com <- red_com
                        new_com[ii] <- 0

                        ii_blue <- match(ii, red_strains)


                        new_com[blue_perm[ii_blue]] <- 1

                        new_com
                    },
                    red_com = red_com, blue_perm = blue_perm,
                    red_strains = red_strains
                    )

                    coms <- apply(coms, 2, paste0, collapse = "")
                    # coms
                    Coms <- c(Coms, coms)
                }
                # Coms

                # Select data
                dat <- Dat %>%
                    filter(community %in% Coms)

                d <- dat %>%
                    mutate(com2 = community) %>%
                    separate(com2, into = paste0("S", 0:8), sep = "") %>%
                    select(-S0)


                cat(paste0(c(red_strains, blue_perm), collapse = ""), "\n")
                Strains <- tibble()
                Ints <- tibble()
                for (i in 1:8) {
                    s <- paste0("S", i)
                    Strains <- Strains %>%
                        bind_rows(tibble(
                            expid = paste0(c(red_strains, blue_perm), collapse = ""),
                            strain = s,
                            effect = d$absorbance[d[[s]] == 1] - d$absorbance[d[[s]] == 0]
                        ))
                    for (j in 2:8) {
                        s2 <- paste("S", j, sep = "")
                        R1 <- d$absorbance[d[[s]] == 1 & d[[s2]] == 0] - d$absorbance[d[[s]] == 0 & d[[s2]] == 0]
                        R2 <- d$absorbance[d[[s]] == 0 & d[[s2]] == 1] - d$absorbance[d[[s]] == 0 & d[[s2]] == 0]

                        if (any(d[[s]] == 1 & d[[s2]] == 1)) {
                            Rint <- d$absorbance[d[[s]] == 1 & d[[s2]] == 1] - d$absorbance[d[[s]] == 0 & d[[s2]] == 0] - R1 - R2
                        } else {
                            Rint <- NA
                        }
                        Ints <- Ints %>%
                            bind_rows(tibble(
                                expid = paste0(c(red_strains, blue_perm), collapse = ""),
                                comb = paste0(sort(c(i, j)), collapse = ""),
                                strain1 = i,
                                strain2 = j,
                                effect = Rint
                            ))
                    }
                }

                Ints
            })
    })



p1 <- Crns %>%
    ggplot(aes(x = effect)) +
    facet_wrap(~comb, scales = "free") +
    geom_density() +
    theme_classic()
p1
ggsave("int_effects.png", p1, width = 6, height = 4)





Crns <- read_tsv("strain_effects.tsv")


p1 <- Crns %>%
    ggplot(aes(x = effect)) +
    facet_wrap(~strain, scales = "free") +
    geom_density() +
    theme_classic()
p1
ggsave("strain_effects.png", p1, width = 6, height = 4)





















