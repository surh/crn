library(tidyverse)
library(coda)

old_opts = options(box.path = "/Users/sur/lab/src/crn/")
box::use(fit_growth/functions)

args <- list()
args$indir <- "/Users/sur/lab/exp/2025/today/single_strain_all_logistic/"
args$outdir <- "/Users/sur/lab/exp/2025/today/"
args$pattern <- "all_logistic"

infiles <- list.files(args$indir, pattern = paste0("_",args$pattern, "_fit.tsv"), 
    full.names = TRUE)

Res <- tibble()
for(f in infiles){
    # f <- infiles[1]
    strain <- basename(f) %>% str_remove(paste0("_",args$pattern, "_fit.tsv"))
    cat(strain, "\n")

    Dat <- read_tsv(f)

    Dat <- Dat %>%
        mutate(r_28 = exp(r_0 + `b_rt[28]`)) %>%
        mutate(r_32 = exp(r_0 + `b_rt[32]`)) %>%
        mutate(K_28 = exp(K_0 + `b_kt[28]`)) %>%
        mutate(K_32 = exp(K_0 + `b_kt[32]`)) %>%
        mutate(delta_r = r_32 - r_28) %>%
        mutate(delta_K = K_32 - K_28) %>%
        select(iteration, chain, r_28, r_32, delta_r, K_28, K_32, delta_K)

    # Dat %>%
    #     ggplot(aes(x = iteration, y = r_28)) +
    #     geom_line(aes(col = factor(chain))) +
    #     theme_classic()

    # Dat %>%
    #     ggplot(aes(x = r_28)) +
    #     geom_histogram(bins = 100) +
    #     # facet_wrap(~chain) +
    #     geom_vline(xintercept = 1.44, linetype = "dashed", color = "red") +
    #     geom_vline(xintercept = 1.03, linetype = "dashed", color = "blue") +
    #     geom_vline(xintercept = 0.598, linetype = "dashed", color = "green") +
    #     theme_classic()
        
    res <- bind_rows(functions$calculate_stats(Dat = Dat, column = "r_28"),
        functions$calculate_stats(Dat = Dat, column = "r_32"),
        functions$calculate_stats(Dat = Dat, column = "delta_r"),
        functions$calculate_stats(Dat = Dat, column = "K_28"),
        functions$calculate_stats(Dat = Dat, column = "K_32"),
        functions$calculate_stats(Dat = Dat, column = "delta_K")) %>%
        mutate(strain = strain) %>%
        select(strain, everything())

    Res <- bind_rows(Res, res)

}
Res %>%
    select(-hpd10, -hpd90) %>%
    print(n = 200)
outfile <- file.path(args$outdir, paste0("single_strain_", args$pattern, "_summary.tsv"))
write_tsv(Res, outfile)
Res <- Res %>%
    filter(Rhat < 1.01) %>%
    print()



# Plotting the results
p1 <- Res %>%
    filter(parameter %in% c("r_28", "r_32")) %>%
    pivot_longer(parameter, names_to = NULL, values_to = "temp") %>%
    ggplot(aes(x = temp, y = mean)) +
    facet_grid(. ~ strain, scales = "free_y") +
    geom_linerange(aes(ymin = hpd05, ymax = hpd95), linewidth = 1, color = "darkgrey") +
    geom_linerange(aes(ymin = hpd10, ymax = hpd90), linewidth = 3, color = "darkgrey") +
    geom_point(size = 8, color = "#dd1c77") +
    geom_point(aes(y = median), size = 6, color = "#c994c7") +
    geom_point(aes(y = mode), size = 4, color = "#7a0177") +
    labs(y = "Growth rate (1/h)", x = "Temperature (ºC)") +
    theme_classic()
p1
outfile <- file.path(args$outdir, paste0("single_strain_", args$pattern, "_r.png"))
ggsave(outfile, p1, width = 7, height = 4)
outfile <- file.path(args$outdir, paste0("single_strain_", args$pattern, "_r.svg"))
ggsave(outfile, p1, width = 7, height = 4)

p1 <- Res %>%
    filter(parameter %in% c("K_28", "K_32")) %>%
    pivot_longer(parameter, names_to = NULL, values_to = "temp") %>%
    ggplot(aes(x = temp, y = mean)) +
    facet_grid(. ~ strain, scales = "free_y") +
    geom_linerange(aes(ymin = hpd05, ymax = hpd95), linewidth = 1, color = "darkgrey") +
    geom_linerange(aes(ymin = hpd10, ymax = hpd90), linewidth = 3, color = "darkgrey") +
    geom_point(size = 8, color = "#dd1c77") +
    geom_point(aes(y = median), size = 6, color = "#c994c7") +
    geom_point(aes(y = mode), size = 4, color = "#7a0177") +
    labs(y = "Carrying capacity (K)", x = "Temperature (ºC)", ) +
    theme_classic() 
p1
outfile <- file.path(args$outdir, paste0("single_strain_", args$pattern, "_K.png"))
ggsave(outfile, p1, width = 7, height = 4)
outfile <- file.path(args$outdir, paste0("single_strain_", args$pattern, "_K.svg"))
ggsave(outfile, p1, width = 7, height = 4)