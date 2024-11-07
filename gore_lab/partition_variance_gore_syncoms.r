# setwd("/Users/sur/lab/exp/2024/today2")
library(tidyverse)
library(lme4)
# library(brms)

#' Code to partition the variance into VPlas and VGen. Code from
#' https://github.com/devillemereuil/CodePartReacNorm from
#' 
#'   Villemereuil, Pierre de, and Luis-Miguel Chevin.
#'   “Partitioning the Phenotypic Variance of Reaction Norms,” September 1,
#'   2023. https://ecoevorxiv.org/repository/view/5891/.

#' Set global parameters
args <- list()
args$input <- "data_from_gore_for_models/day1_data.tsv"


#' Read data and plot
Dat <- read_tsv(args$input,
                col_types = cols(Day = col_number(),
                                 Rep = col_character(),
                                 OD = col_number(), 
                                 S0 = col_number(),
                                 ID = col_character(),
                                 Nutrient = col_number(),
                                 Community = col_character()))
Dat 

p1 <- Dat %>%
  ggplot(aes(x = Nutrient, y = OD,
             group = interaction(Community,Rep, sep = "_")) ) +
  geom_line(aes(col = factor(S0))) +
  scale_color_brewer(type = "seq", palette = "YlOrRd") + 
  AMOR::theme_blackbox() +
  theme(panel.background = element_rect(fill = "grey"))
p1

#' We also create a nutrient squared column for the polynomial models and
#' a factor version of nutrient for the character state model
Dat <- Dat %>%
  mutate(Nutrient_sq = Nutrient * Nutrient,
         Nutrient_f = factor(Nutrient))

#' # Model
#' We are going to try 3 basic polynomial models with degrees 0-2

mp0 <- lmer(OD ~ 1 + (1 | Community), data = Dat)
summary(mp0)

mp1 <- lmer(OD ~ 1 + Nutrient + (1 + Nutrient | Community), data = Dat)
summary(mp1)

mp2 <- lmer(OD ~ 1 + Nutrient + Nutrient_sq + (1 + Nutrient + Nutrient_sq | Community), data = Dat)
summary(mp2)

#' We also fit a character state model
mc <- lmer(OD ~ 0 + Nutrient_f + (0 + Nutrient_f | Community),
           data = Dat)
summary(mc)

#' Compare all models with AIC, BIC
AIC(mp0, mp1, mp2, mc)
BIC(mp0, mp1, mp2, mc)

#' ## Extract params from polynomial models
partition_variance_polynomial <- function(mp,
                                          pheno_name = "OD"){
  #' Design matrix
  design_mat <- model.matrix(mp)
  
  #' Variance-covariance matrix
  varcov <- VarCorr(mp)[["Community"]]
  attr(varcov, "stddev") <- attr(varcov, "correlation") <- NULL
  
  # Average (fixed) effects
  out <- summary(mp)[["coefficients"]]
  coefs <- out[ , "Estimate"]
  se    <- out[ , "Std. Error"]
  
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
  V_Gen <- sum( diag( (1/nrow(design_mat)) * (t(design_mat) %*% design_mat) %*% varcov ) )
  
  
  # V_Tot <- V_Plas + V_Gen + V_Res
  V_Phen <- var(model.frame(mp)[,pheno_name])
  
  Pi <- (coefs^2 * diag(varcov_design) - se^2) / V_Plas
  Gamma <- ( ( ( t(design_mat) %*% design_mat ) / nrow(design_mat) ) * varcov ) / V_Gen
  
  return(list(design_mat = design_mat,
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
              Gamma = Gamma))
}

partition_variance_polynomial(mp = mp0, pheno_name = "OD")
partition_variance_polynomial(mp = mp1, pheno_name = "OD")
partition_variance_polynomial(mp = mp2, pheno_name = "OD")


#' ## Extract params from character state model











