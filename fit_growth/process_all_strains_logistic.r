library(tidyverse)
library(coda)

old_opts = options(box.path = "/Users/sur/lab/src/crn/")
box::use(fit_growth/functions)

args <- list()
args$indir <- "/Users/sur/lab/exp/2025/today3/all_strains_logistic/"
args$outdir <- "/Users/sur/lab/exp/2025/today3/"


Dat <- read_tsv(file.path(args$indir, "all_logistic_fit.tsv"))

#' Identify variables and reshape data for tidyverse
dat <- Dat %>% 
    pivot_longer(-c(iteration, chain, lp, n_steps, is_accept, 
        acceptance_rate, log_density, hamiltonian_energy, 
        hamiltonian_energy_error, max_hamiltonian_energy_error, 
        tree_depth, numerical_error, step_size, nom_step_size), 
        names_to = "parameter", values_to = "value") %>%
        select(iteration, chain, parameter, value) %>%
        mutate(type = parameter,
                strain = NA,
                temp = NA) %>%
        mutate(type = replace(type, str_detect(parameter,  "^b_rs\\["), "b_rs")) %>%
        mutate(type = replace(type, str_detect(parameter,  "^b_rt\\["), "b_rt")) %>%
        mutate(type = replace(type, str_detect(parameter,  "^b_rst\\["), "b_rst")) %>%
        mutate(type = replace(type, str_detect(parameter,  "^b_rb\\["), "b_rb")) %>%
        mutate(type = replace(type, str_detect(parameter,  "^b_ks\\["), "b_ks")) %>%
        mutate(type = replace(type, str_detect(parameter,  "^b_kt\\["), "b_kt")) %>%
        mutate(type = replace(type, str_detect(parameter,  "^b_kst\\["), "b_kst")) %>%
        mutate(type = replace(type, str_detect(parameter,  "^b_kb\\["), "b_kb"))
dat$strain[ dat$type == "b_rs" ] <- dat$parameter[ dat$type == "b_rs" ] %>%
    str_remove("b_rs\\[String7\\(\"") %>%
    str_remove("\"\\)\\]$")
dat$strain[ dat$type == "b_ks" ] <- dat$parameter[ dat$type == "b_ks" ] %>%
    str_remove("b_ks\\[String7\\(\"") %>%
    str_remove("\"\\)\\]$")

dat$strain[ dat$type == "b_rst" ] <- dat$parameter[ dat$type == "b_rst" ] %>%
    str_remove("b_rst\\[\"") %>%
    str_remove("_[\\d]{2}\\.0\"\\]$")
dat$strain[ dat$type == "b_kst" ] <- dat$parameter[ dat$type == "b_kst" ] %>%
    str_remove("b_kst\\[\"") %>%
    str_remove("_[\\d]{2}\\.0\"\\]$")

dat$temp[ dat$type == "b_rt" ] <- dat$parameter[ dat$type == "b_rt" ] %>%
    str_remove("b_rt\\[") %>%
    str_remove("\\]$")
dat$temp[ dat$type == "b_kt" ] <- dat$parameter[ dat$type == "b_kt" ] %>%
    str_remove("b_kt\\[") %>%
    str_remove("\\]$")

dat$temp[ dat$type == "b_rst" ] <- dat$parameter[ dat$type == "b_rst" ] %>%
    str_remove("b_rst\\[\"ST[\\d]{5}_") %>%
    str_remove("\"\\]$")
dat$temp[ dat$type == "b_kst" ] <- dat$parameter[ dat$type == "b_kst" ] %>%
    str_remove("b_kst\\[\"ST[\\d]{5}_") %>%
    str_remove("\"\\]$")

# table(dat$strain, useNA = "always")
# table(dat$temp, useNA = "always")
# table(dat$strain,dat$temp, useNA = "always")

#' Reformat data to get r and K for each strain and temperature
dat <- dat %>%
    # filter(iteration < 600) %>%
    # filter(chain == 1) %>%
    filter(type %in% c("r_0", "b_rs", "b_rt", "b_rst", 
        "K_0", "b_ks", "b_kt", "b_kst")) %>%
    # print(n = 1000)
    group_split(.$iteration, .$chain)  %>%
    map_dfr(function(d){
        strains <- unique(d$strain[!is.na(d$strain)])
        temps <- unique(d$temp[!is.na(d$temp)])
        iteration <- unique(d$iteration)
        chain <- unique(d$chain)
        # print(strains)
        # print(d, n = 10000)

        res <- tibble()
        for(strain in strains){
            for(temp in temps){
                r <- d$value[d$type == "r_0"] + 
                    d$value[d$type == "b_rs" & strain == d$strain] +
                    d$value[d$type == "b_rt" & temp == d$temp] +
                    d$value[d$type == "b_rst" & strain == d$strain & temp == d$temp]

                K <- d$value[d$type == "K_0"] + 
                    d$value[d$type == "b_ks" & strain == d$strain] +
                    d$value[d$type == "b_kt" & temp == d$temp] +
                    d$value[d$type == "b_kst" & strain == d$strain & temp == d$temp]

                res <- bind_rows(res,
                    tibble(strain = strain,
                            temp = temp,
                            r = r,
                            K = K))

            }
        }
        res$chain <- chain
        res$iteration <- iteration
        res
    })
dat

#' Calculate statistics for growth rate and carrying capacity for each strain
#' and temperature
strains <- unique(dat$strain)
temps <- unique(dat$temp)
Res <- tibble()
for(st in strains){
    for(tp in temps){
        Res <- bind_rows(Res,    
            functions$calculate_stats(Dat = dat %>%
                filter(strain == st, temp == tp),
                column = "r") %>%
                mutate(strain = st, temp = tp),
            functions$calculate_stats(Dat = dat %>%
                filter(st == strain, temp == tp),
                column = "K") %>%
                mutate(strain = st, temp = tp))

    }
}
Res <- Res %>%
    select(strain, temp, Rhat, everything())
Res
write_tsv(Res, file.path(args$outdir, "all_strains_logistic_summary.tsv"))

Res <- Res %>%
    filter(Rhat < 1.01) %>%
    print()

# Res %>% filter(parameter == "r_28") %>% print(n=100)
# Res %>% filter(parameter == "K") %>% print(n=100)

p1 <- Res %>%
    filter(parameter %in% c("r")) %>%
    # pivot_longer(parameter, names_to = NULL, values_to = "temp") %>%
    ggplot(aes(x = temp, y = mean)) +
    facet_grid(. ~ strain, scales = "free_y") +
    geom_linerange(aes(ymin = hpd05, ymax = hpd95), linewidth = 1, color = "darkgrey") +
    geom_linerange(aes(ymin = hpd10, ymax = hpd90), linewidth = 3, color = "darkgrey") +
    geom_point(size = 6, color = "#dd1c77") +
    geom_point(aes(y = median), size = 5, color = "#c994c7") +
    labs(y = "Growth rate (1/h)", x = "Temperature (ºC)") +
    theme_classic()
p1
outfile <- file.path(args$outdir, "all_strains_logistic_r.png")
ggsave(outfile, p1, width = 7, height = 4)
outfile <- file.path(args$outdir, "all_strains_logistic_r.svg")
ggsave(outfile, p1, width = 7, height = 4)


#! Need to adapt to different K per temperature!!!
p1 <- Res %>%
    filter(parameter %in% c("K")) %>%
    # pivot_longer(parameter, names_to = NULL, values_to = "temp") %>%
    ggplot(aes(x = temp, y = mean)) +
    facet_grid(. ~ strain, scales = "free_y") +
    geom_linerange(aes(ymin = hpd05, ymax = hpd95), linewidth = 1, color = "darkgrey") +
    geom_linerange(aes(ymin = hpd10, ymax = hpd90), linewidth = 3, color = "darkgrey") +
    geom_point(size = 6, color = "#dd1c77") +
    geom_point(aes(y = median), size = 5, color = "#c994c7") +
    labs(y = "Growth rate (1/h)", x = "Temperature (ºC)") +
    theme_classic()
p1
outfile <- file.path(args$outdir, "all_strains_logistic_K.png")
ggsave(outfile, p1, width = 7, height = 4)
outfile <- file.path(args$outdir, "all_strains_logistic_K.svg")
ggsave(outfile, p1, width = 7, height = 4)