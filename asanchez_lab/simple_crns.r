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




# There are 70 combinations of 4 species
# For each combination the plan is compare sampling random 16 communities
# with the crn approach and the dropout approach. Should include the exhaustive
# apporach as well. I don't have the interaction estimates from the paper,
# so the exhaustive approach will make do as an alternative

# First sample random 16 communities 100 times
n_coms <- 16
set.seed(1503)
Res <- NULL
for (i in 1:1680){
    cat(i, "\n")
    coms <- unique(Dat$community)
    dat <- Dat %>%
        filter(community %in% sample(coms, size = n_coms, replace = FALSE))
    # dat

    # Model
    # We are going to try 3 basic polynomial models with degrees 0-2
    mp0 <- lmer(absorbance ~ 1 + (1 | community), data = dat)
    # summary(mp0)

    mp1 <- lmer(absorbance ~ 1 + dilution_factor + (1 + dilution_factor | community), data = dat)
    # summary(mp1)

    mp2 <- lmer(absorbance ~ 1 + dilution_factor + dilution_factor_sq + (1 + dilution_factor + dilution_factor_sq | community), data = dat)
    # summary(mp2)
    
    # We also fit a character state model. DON't FIT, can't extract
    # params yet
    # mc <- lmer(absorbance ~ 0 + dilution_factor_f + (0 + dilution_factor_f | community),
    #        data = dat)
    # summary(mc)

    # Select the best model
    m_ii <- which.min(AIC(mp0, mp1, mp2)$AIC)
    m_selected <- list(mp0, mp1, mp2)[[m_ii]]
    # BIC(mp0, mp1, mp2, mc)


    res <- partition_variance_polynomial(
        mp = m_selected,
        pheno_name = "absorbance",
        com_name = "community"
    )

    Res <- Res %>%
        bind_rows( tibble(
        expid = paste0("r",i),
        m_type = c("mp0", "mp1", "mp2")[m_ii],
        V_Phen = res$V_Phen,
        V_Tot = res$V_Tot,
        V_Plas = res$V_Plas,
        V_Gen = res$V_Gen,
        V_Res = res$V_Res
    ))

    rm(mp0, mp1, mp2, m_selected)
}
Res
Res %>%
    write_tsv(file.path(args$outdir, "random_draws_crns.tsv"))



# Full data
# Model
# We are going to try 3 basic polynomial models with degrees 0-2
mp0 <- lmer(absorbance ~ 1 + (1 | community), data = Dat)
# summary(mp0)

mp1 <- lmer(absorbance ~ 1 + dilution_factor + (1 + dilution_factor | community), data = Dat)
# summary(mp1)

mp2 <- lmer(absorbance ~ 1 + dilution_factor + dilution_factor_sq + (1 + dilution_factor + dilution_factor_sq | community), data = Dat)
# summary(mp2)


AIC(mp0, mp1, mp2)

# Select the best model
m_ii <- which.min(AIC(mp0, mp1, mp2)$AIC)
m_selected <- list(mp0, mp1, mp2)[[m_ii]]


res <- partition_variance_polynomial(
    mp = m_selected,
    pheno_name = "absorbance",
    com_name = "community"
)

tibble(
    m_type = c("mp0", "mp1", "mp2")[m_ii],
    V_Phen = res$V_Phen,
    V_Tot = res$V_Tot,
    V_Plas = res$V_Plas,
    V_Gen = res$V_Gen,
    V_Res = res$V_Res
) %>%
    print() %>%
    write_tsv(file.path(args$outdir, "full_data_crns.tsv"))
summary(Res)



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

                # Model
                # We are going to try 3 basic polynomial models with degrees 0-2
                mp0 <- lmer(absorbance ~ 1 + (1 | community), data = dat)
                # summary(mp0)

                mp1 <- lmer(absorbance ~ 1 + dilution_factor + (1 + dilution_factor | community), data = dat)
                # summary(mp1) dilution_factor | community), data = dat, control = lmerControl(check.nobs.vs.nRE = "ignore"))
                # summary(mp1)

                mp2 <- lmer(absorbance ~ 1 + dilution_factor + dilution_factor_sq + (1 + dilution_factor + dilution_factor_sq | community), data = dat)
                # summary(mp2)

                # We also fit a character state model. DON't FIT, can't extract
                # params yet
                # mc <- lmer(absorbance ~ 0 + dilution_factor_f + (0 + dilution_factor_f | community),
                #        data = dat)
                # summary(mc)

                # Select the best model
                m_ii <- which.min(AIC(mp0, mp1, mp2)$AIC)
                m_selected <- list(mp0, mp1, mp2)[[m_ii]]
                # BIC(mp0, mp1, mp2, mc)


                res <- partition_variance_polynomial(
                    mp = m_selected,
                    pheno_name = "absorbance",
                    com_name = "community"
                )


                Res <- tibble(
                    expid = paste0(c(red_strains, blue_perm), collapse = ""),
                    m_type = c("mp0", "mp1", "mp2")[m_ii],
                    V_Phen = res$V_Phen,
                    V_Tot = res$V_Tot,
                    V_Plas = res$V_Plas,
                    V_Gen = res$V_Gen,
                    V_Res = res$V_Res
                )

                rm(mp0, mp1, mp2, m_selected)
                Res
            })

        Res.per.red
    })

Crns
Crns %>%
    write_tsv(file.path(args$outdir, "f2like_crns.tsv"))






















