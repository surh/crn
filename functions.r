#' ## Extract params from polynomial models
partition_variance_polynomial <- function(mp,
                                          pheno_name = "OD", 
                                          com_name = "Community"){

    #' Design matrix
    design_mat <- model.matrix(mp)

    #' Variance-covariance matrix
    varcov <- VarCorr(mp)[[com_name]]
    attr(varcov, "stddev") <- attr(varcov, "correlation") <- NULL

    # Average (fixed) effects
    out <- summary(mp)[["coefficients"]]
    coefs <- out[, "Estimate"]
    se <- out[, "Std. Error"]

    # Residual variance
    V_Res <- attr(VarCorr(mp), "sc")^2

    # Variance-covariance of environment
    varcov_design <- cov(design_mat)

    # For a polynomial model we have an unbiased estimator for V_plas
    # \hat{V}_{Plas} = \bar\theta^TX\theta - Tr(S_\theta X)
    V_Plas <- coefs %*% varcov_design %*% coefs - se %*% varcov_design %*% se
    V_Plas <- as.numeric(V_Plas)

    #' For V_gen we have:
    #' V_{Gen} = E_\epsilon(x^T\Theta x) = \bar{x}^T\Theta\bar{x} + Tr(\Theta X)
    #' Is this equivalent?
    V_Gen <- sum(diag((1 / nrow(design_mat)) * (t(design_mat) %*% design_mat) %*% varcov))


    # V_Tot <- V_Plas + V_Gen + V_Res
    V_Phen <- var(model.frame(mp)[, pheno_name])

    Pi <- (coefs^2 * diag(varcov_design) - se^2) / V_Plas
    Gamma <- (((t(design_mat) %*% design_mat) / nrow(design_mat)) * varcov) / V_Gen

    return(list(
        design_mat = design_mat,
        varcov_design = varcov_design,
        varcov = varcov,
        coefs = coefs,
        se = se,
        V_Phen = V_Phen,
        V_Tot = V_Plas + V_Gen + V_Res,
        V_Plas = V_Plas,
        V_Gen = V_Gen,
        V_Res = V_Res,
        Pi = Pi,
        Gamma = Gamma
    ))
}
