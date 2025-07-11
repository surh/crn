library(tidyverse)
library(coda)

old_opts = options(box.path = "/Users/sur/lab/src/crn/")
box::use(fit_growth/functions)

args <- list()
args$indir <- "/Users/sur/lab/exp/2025/today3/single_strain_all_logistic/"
args$outdir <- "/Users/sur/lab/exp/2025/today3/"

infiles <- list.files(args$indir, pattern = "all_logistic_fit.tsv", full.names = TRUE)

Res <- tibble()
for(f in infiles){
    # f <- infiles[1]
    strain <- basename(f) %>% str_remove("_all_logistic_fit.tsv")
    cat(strain, "\n")

    Dat <- read_tsv(f)


    Dat <- Dat %>%
        mutate(r_28 = r_0 + `b_rt[28]`) %>%
        mutate(r_32 = r_0 + `b_rt[32]`) %>%
        mutate(delta_r = r_32 - r_28) %>%
        select(iteration, chain, K, r_28, r_32, delta_r)
        
    res <- bind_rows(functions$calculate_stats(Dat = Dat, column = "r_28"),
        functions$calculate_stats(Dat = Dat, column = "r_32"),
        functions$calculate_stats(Dat = Dat, column = "delta_r"),
        functions$calculate_stats(Dat = Dat, column = "K")) %>%
        mutate(strain = strain) %>%
        select(strain, everything())

    Res <- bind_rows(Res, res)

}
Res
Res <- Res %>%
    filter(Rhat < 1.01) %>%
    print()

# Res %>% filter(parameter == "r_28") %>% print(n=100)
# Res %>% filter(parameter == "K") %>% print(n=100)

p1 <- Res %>%
    filter(parameter %in% c("r_28", "r_32")) %>%
    pivot_longer(parameter, names_to = NULL, values_to = "temp") %>%
    ggplot(aes(x = temp, y = mean)) +
    facet_grid(. ~ strain, scales = "free_y") +
    geom_linerange(aes(ymin = hpd05, ymax = hpd95), linewidth = 1, color = "darkgrey") +
    geom_linerange(aes(ymin = hpd10, ymax = hpd90), linewidth = 3, color = "darkgrey") +
    geom_point(size = 6, color = "#dd1c77") +
    geom_point(aes(y = median), size = 5, color = "#c994c7") +
    labs(y = "Growth rate (1/h)", x = "Temperature (ºC)") +
    theme_classic()
p1
outfile <- file.path(args$outdir, "single_strain_all_logistic_r.png")
ggsave(outfile, p1, width = 7, height = 4)
outfile <- file.path(args$outdir, "single_strain_all_logistic_r.svg")
ggsave(outfile, p1, width = 7, height = 4)


#! Need to adapt to different K per temperature!!!
p1 <- Res %>%
    filter(parameter == "K") %>%
    ggplot(aes(x = temp, y = mean)) +
    facet_grid(. ~ strain, scales = "free_y") +
    geom_linerange(aes(ymin = hpd05, ymax = hpd95), linewidth = 1, color = "darkgrey") +
    geom_linerange(aes(ymin = hpd10, ymax = hpd90), linewidth = 3, color = "darkgrey") +
    geom_point(size = 6, color = "#dd1c77") +
    geom_point(aes(y = median), size = 5, color = "#c994c7") +
    labs(y = "Carrying capacity (K)", x = "Temperature (ºC)") +
    theme_classic()
p1
outfile <- file.path(args$outdir, "single_strain_bytemp_logistic_K.png")
ggsave(outfile, p1, width = 7, height = 4)
outfile <- file.path(args$outdir, "single_strain_bytemp_logistic_K.svg")
ggsave(outfile, p1, width = 7, height = 4)