conflicted::conflict_prefer("pack",     "tidyr")
conflicted::conflict_prefer("unpack",   "tidyr")
conflicted::conflict_prefer("expand",   "tidyr")
conflicted::conflict_prefer("extract",  "tidyr")
conflicted::conflict_prefer("filter",   "dplyr")
conflicted::conflict_prefer("lag",      "dplyr")
conflicted::conflict_prefer("chol2inv", "Matrix")

library(tidyverse)
# library(mvtnorm)
library(lme4)
# library(furrr)
# library(patchwork)
# library(here)

## Function to run the reaction norm model on the data
# Args: - df: A dataset containing the phenotypic data
# Value: The model fit for the data
fit_lme_rn <- function(df) {
  # Formatting the data
  df[["Env_Sq"]] <- df[["Env"]]^2
  
  # Fitting the model
  lmer(Phen ~ Env + Env_Sq + (1 + Env + Env_Sq | Genotype),
       data = df)
}

## Function to run the character-state model on the data
# Args: - df: A dataset containing the phenotypic data
# Value: The model fit for the data
fit_lme_cs <- function(df) {
  # Formatting the data
  df[["Env_fac"]] <- as_factor(round(df[["Env"]], digits = 1))
  
  # Fitting the model
  lmer(Phen ~ 0 + Env_fac + (0 + Env_fac | Genotype),
       data = df)
}

## Function to extract the parameters of interest from the RN lme4 fit
# Args: - model: A lme fit yielded by fit_lme_rn
# Values: List of relevant parameters
extract_params_rn <- function(model) {
  # Variance-covariance matrix
  mat_vcv <- VarCorr(model)[["Genotype"]]
  attr(mat_vcv, "stddev") <- attr(mat_vcv, "correlation") <- NULL
  
  # Average (fixed) effects
  out   <- summary(model)[["coefficients"]]
  coefs <- out[ , "Estimate"]
  se    <- out[ , "Std. Error"]
  
  # Residual variance
  res_var <- attr(VarCorr(model), "sc")^2
  
  # Returning the list of parameters
  tibble(a    = coefs[1],
         b    = coefs[2],
         c    = coefs[3],
         a_se = se[1],
         b_se = se[2],
         c_se = se[3],
         Va   = mat_vcv[1, 1],
         Vb   = mat_vcv[2, 2],
         Vc   = mat_vcv[3, 3],
         Cab  = mat_vcv[1, 2],
         Cac  = mat_vcv[1, 3],
         Cbc  = mat_vcv[2, 3],
         Vr   = res_var)
}

## Function to extract the parameters of interest from the CS lme4 fit
# Args: - model: A lme fit yielded by fit_lme_cs
# Values: List of relevant parameters
extract_params_cs <- function(model) {
  # Genotypic VCV
  G <- VarCorr(model)[["Genotype"]]
  attr(G, "stddev") <- NULL
  attr(G, "correlation") <- NULL
  
  # Average (fixed) effects
  out   <- summary(model)[["coefficients"]]
  var_p <- var(predict(model, re.form = ~ 0)) - mean(out[, "Std. Error"]^2)
  
  # Residual variance
  res_var <- attr(VarCorr(model), "sc")^2
  
  # Returning the list of parameters
  tibble(V_Plas_CS   = var_p,
         V_Gen_CS    = mean(diag(G)),
         V_Res_CS    = res_var)
}

## Function to compute the expected genotypic variance
# Args: - vcv: the variance-covariance matrix between the RN parameters
#       - mat_env: the "design matrix" of the reaction norm
# Values: The expected genotypic variance
compute_var_geno <- function(vcv, mat_env) {
  # Sum of the element-wise multiplication of products with mat_env
  # with the elements of the variance-covariance matrix vcv
  # Equivalent (but faster) than averaging over the t(E[i, ]) %*% vcv %*% E[i, ] products
  sum(diag((1/nrow(mat_env)) * (t(mat_env) %*% mat_env) %*% vcv))
}

env <- seq(-1, 1, length.out = 3)

args <- list()
args$input <- "data_from_gore_for_models/day1_data.tsv"
Dat <- read_tsv(args$input,
                col_types = cols(Day = col_number(),
                                 Rep = col_character(),
                                 OD = col_number(), 
                                 S0 = col_number(),
                                 ID = col_character(),
                                 Nutrient = col_number(),
                                 Community = col_character()))
Dat <- Dat %>%
  rename(Phen = OD, Genotype = Community, Env = Nutrient)
Dat



mp2.geno <- fit_lme_rn(Dat)
mc.geno <- fit_lme_cs(Dat)

## Saving the design matrix
mat_env <- model.matrix(mp2.geno)
## Compute the variance-covariance matrix of the environmental values
vcv_env <- cov(mat_env)

Res <- tibble(Day = 1, Data = list(tibble(Dat))) %>%
  mutate(Params_RN = extract_params_rn(mp2.geno)) %>%
  tidyr::unpack(Params_RN) %>%
  mutate(Params_CS = extract_params_cs(mc.geno)) %>%
  tidyr::unpack(Params_CS)
Res

Res |>
  rowwise() |>
  mutate(V_Plas_RN    = as.vector(t(c(a, b, c)) %*% vcv_env %*% c(a, b, c)) -
           as.vector(t(c(a_se, b_se, c_se)) %*% vcv_env %*% c(a_se, b_se, c_se)),
         V_Gen_RN     = compute_var_geno(matrix(c(Va, Cab, Cac, Cab, Vb, Cbc, Cac, Cbc, Vc),
                                                ncol = 3, nrow = 3),
                                         mat_env),
         Pi_RN_b  = ((b^2 * vcv_env[2, 2] - b_se^2) / V_Plas_RN),
         Pi_RN_c  = ((c^2 * vcv_env[3, 3] - c_se^2) / V_Plas_RN),
         Gamma_RN_a   = Va / V_Gen_RN,
         Gamma_RN_b   = Vb * mean(env^2) / V_Gen_RN,
         Gamma_RN_c   = Vc * mean(env^4) / V_Gen_RN,
         Gamma_RN_ac  = 2 * Cac * mean(env^2) / V_Gen_RN,
         V_Res_RN     = Vr,
         V_Tot_RN     = V_Plas_RN + V_Gen_RN + V_Res_RN,
         V_Tot_CS     = V_Plas_CS + V_Gen_CS + V_Res_CS,
         V_Phen       = var(Data[["Phen"]])) |>
  ungroup() %>%
  select(Day, contains("V_"), starts_with("Pi"), starts_with("Gamma")) %>%
  print(width = 1e4)



