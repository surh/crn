#' ## Extract params from polynomial models
partition_variance_polynomial <- function(mp,
                                          pheno_name = "AUC", 
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

read_single_experiment <- function(od600_file, timepoints_file, syncoms_file = NULL,
                                   mpn_file = NULL, type = "syncom",
                                   batch_name =  basename(dirname(od600_file))){
  
  # od600_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/od600.tsv"
  # timepoints_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/timepoints.tsv"
  # syncoms_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/syncoms.tsv"
  # mpn_file <- "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/mpn.tsv"
  # type <- "syncom"
  # batch_name <- "NS1"
  
  
  # od600_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/NS1/od600.tsv"
  # timepoints_file = "/Users/sur/lab/data/2024_rhizo_pilot_syncom_NS/single_strains/NS1/timepoints.tsv"
  # type <- "strain"
  # batch_name <- "SS_NS1"

  if(type == "syncom"){
    id_cols <- c("Community", "temp")
  }else if(type == "strain"){
    id_cols <- c("strain", "temp")
  }
  
  
  #' Read and proces OD
  OD <- read_tsv(od600_file)
  Timepoints <- read_tsv(timepoints_file)
  
  Dat <- OD %>%
    pivot_longer(-id_cols,
                 names_to = "timepoint",
                 values_to = "OD600") %>%
    left_join(Timepoints, by = c("timepoint", "temp")) %>%
    mutate(OD600 = (OD600 - blank_od) * od_factor) %>%
    mutate(OD600 = replace(OD600, timepoint == "t_0", 0.001)) 
  
  
  #' Add MPN if available
  if(!is.null(mpn_file) && file.exists(mpn_file)){
    MPN <- read_tsv(mpn_file, na = c("", "NA", "Error"))
    
    Dat <- MPN %>%
      pivot_longer(-id_cols,
                   names_to = "timepoint",
                   values_to = "MPN") %>%
      right_join(Dat, by = c(id_cols, "timepoint")) 
  }
  
  if(!is.null(batch_name)){
    Dat[["batch"]] <- batch_name
  }
  
  
  return(Dat)
  
}