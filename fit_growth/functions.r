box::use(magrittr[...])


#' @export
calculate_stats <- function(Dat, column = "r"){

    par.mcmc <- coda::as.mcmc(Dat[[column]])
    hpd <- coda::HPDinterval(par.mcmc, prob = 0.9)
    hpd05 <- hpd[1]
    hpd95 <- hpd[2]
    hpd <- coda::HPDinterval(par.mcmc, prob = 0.8)
    hpd10 <- hpd[1]
    hpd90 <- hpd[2]


    Rhat <- Dat %>%
        dplyr::select(!!column, chain, iteration) %>%
        tidyr::pivot_wider(names_from = chain, values_from = !!column) %>%
        dplyr::select(-iteration) %>%
        as.matrix() %>%
        rstan::Rhat()

    dplyr::tibble(parameter = column,
           mean = mean(as.numeric(Dat[[column]]), na.rm = TRUE),
           sd = stats::sd(as.numeric(Dat[[column]]), na.rm = TRUE),
           median = stats::median(as.numeric(Dat[[column]]), na.rm = TRUE),
           hpd05 = hpd05,
           hpd10 = hpd10,
           hpd90 = hpd90,
           hpd95 = hpd95,
           Rhat = Rhat,
           q05 = stats::quantile(as.numeric(Dat[[column]]), 0.05, na.rm = TRUE),
           q10 = stats::quantile(as.numeric(Dat[[column]]), 0.1, na.rm = TRUE),
           q90 = stats::quantile(as.numeric(Dat[[column]]), 0.9, na.rm = TRUE),
           q95 = stats::quantile(as.numeric(Dat[[column]]), 0.95, na.rm = TRUE)

    )
}


