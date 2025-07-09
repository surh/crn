library(tidyverse)
library(coda)

old_opts = options(box.path = "/Users/sur/lab/src/crn/")
box::use(fit_growth/functions)

args <- list()
args$indir <- "/Users/sur/lab/exp/2025/today3/single_strain_bytemp_logistic/"
args$outdir <- "/Users/sur/lab/exp/2025/today3/"

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
outfile <- file.path(args$outdir, "single_strain_bytemp_logistic_r.png")
ggsave(outfile, p1, width = 6, height = 4)
outfile <- file.path(args$outdir, "single_strain_bytemp_logistic_r.svg")
ggsave(outfile, p1, width = 6, height = 4)


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