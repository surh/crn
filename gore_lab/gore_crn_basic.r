# setwd("/Users/sur/lab/exp/2024/today2")
library(tidyverse)
library(lme4)
library(brms)


#' Reference:
#' https://github.com/devillemereuil/TutoPartReacNorm/blob/main/tuto_discrete_environment.R

#' Read data. Clean ID names to remove spaces
Pheno <- read_csv("Matrix_OD.csv")
Pheno

Counts <- read_csv("Matrix_abundance.csv")
names(Counts) <- c("Species", names(Counts)[-1] %>% str_remove(" "))
Counts


Meta <- read_csv("Matrix_reference.csv")
Meta <- Meta %>%
  mutate(ID = str_remove(ID, " "))
Meta

#' Match phenotype data with metadata (i.e. add IDs)
Pheno <- Pheno %>%
  full_join(Meta, by = c("Day", "Rep", "Community", "S0", "Nutrient"))
Pheno

#' Check that no ID appears more than once
if (!all(table(Pheno$ID[!is.na(Pheno$ID)]) == 1)){
  stop("ERROR. At least one ID is repeated in Pheno table")
}

#' Check that all IDs in Meta show up in Pheno
if(length(setdiff(Meta$ID, Pheno$ID[!is.na(Pheno$ID)])) != 0){
  stop("ERROR. At least one 1 from Meta missing in Pheno")
}

#' Generate M IDs for Pheno measurements without IDs
Pheno$ID[is.na(Pheno$ID)] <- paste0("M", 
                                    (max(Pheno$ID[!is.na(Pheno$ID)] %>%
                                           str_remove("^M") %>%
                                           as.numeric()) + 1):nrow(Pheno))
Pheno

#' Recode nutrient as -1 (LN), 0 (MN), and 1 (HN)
Pheno <- Pheno %>%
  rename(Nutrient_orig = Nutrient) %>%
  mutate(Nutrient = NA) %>%
  mutate(Nutrient = replace(Nutrient, Nutrient_orig == "LN", -1)) %>%
  mutate(Nutrient = replace(Nutrient, Nutrient_orig == "MN", 0)) %>%
  mutate(Nutrient = replace(Nutrient, Nutrient_orig == "HN", 1)) %>%
  select(-Nutrient_orig)
Pheno


#' Temporary fix. Create new community with combination of Community and S0
Pheno <- Pheno %>%
  mutate(Com = interaction(Community, S0, sep = "_") %>% as.character()) %>%
  select(-Community) %>%
  rename(Community = Com)
Pheno

#' Model day1 
Dat <- Pheno %>%
  filter(Day == 1)
p1 <- Dat %>%
  ggplot(aes(x = Nutrient, y = OD,
             group = interaction(Community,Rep, sep = "_")) ) +
  geom_line(aes(col = factor(S0))) +
  scale_color_brewer(type = "seq", palette = "YlOrRd") + 
  AMOR::theme_blackbox() +
  theme(panel.background = element_rect(fill = "grey"))
p1

#' * Day 1 reference model without nutrient
m0 <- lmer(OD ~ 1 +  (1  | Community), data = Dat)
summary(m0)

#' * Day 1 simplest (slope and intercept mixed model, classic reaction norm)
m1 <- lmer(OD ~ 1 + Nutrient + (1 + Nutrient | Community), data = Dat)
summary(m1)

#' * Day 1 with S_0
m2 <- lmer(OD ~ 1 + Nutrient + S0 + (1 + Nutrient | Community), data = Dat)
summary(m2)

#' * Day 1 full model Nutrient + S_0 + Nutrient x S_0
m3 <- lmer(OD ~ 1 + Nutrient * S0 + (1 + Nutrient | Community), data = Dat)
summary(m3)


#' Compare all modelos with AIC, BIC
AIC(m0, m1, m2, m3)
BIC(m0, m1, m2, m3)

#' In this particular case it seems like nutrient is important (as expected)
#' and thus models that include this variable perform better. It seems like 
#' S_0 has no measurable effect and thus both BIC and AIC agree to drop it.
#' This in interesting because the figure does show some pattern related to
#' high S_0


#' **To do**: instead of modeling OD directly, we should model \delta OD,
#' meaning the change in OD either from day 0, or from the previous day.


m1.sum <- summary(m1)
m1.sum

#' Re-fit model with brms to get full bayesian benefits

m1.brms <- brm(OD ~ 1 + Nutrient + (1 + Nutrient | Community), data = Dat,
               chains = 4, cores = 4, iter = 2000, thin = 1, warmup = 500)
m1.brms

# Dat %>%
#   mutate(Predict = predict(m1.brms,
#                            re_formula = NA) %>%
#                       as_tibble()) %>%
#   tidyr::unpack(Predict) %>%
#   select(Nutrient,
#          Predict = Estimate,
#          Predict_Low = Q2.5,
#          Predict_Up  = Q97.5) %>%
#   summarise(across(starts_with("Predict"), mean),
#             .by = Nutrient)


# Computing the environmental variance-covariance matrix
design_mat <- model.matrix(OD ~ Nutrient, data = Dat)
cov_X <- cov(design_mat)

# Computing the average of powers of X
M <- (1 / nrow(design_mat)) * t(design_mat) %*% design_mat

## Extracting the estimates from the model ----
beta_hat <- fixef(m1.brms, summary = FALSE)
vcov_hat <- vcov(m1.brms)
G_hat <- VarCorr(m1.brms, summary = FALSE)[["Community"]][["cov"]] %>%
  apply(1, \(mat) { mat }, simplify = FALSE)

# Setting up the correcting variance accounting for the uncertainty issue (Eq. 43)
var_uncert <- sum(cov_X * vcov_hat)
var_uncert

V_I_hat  <- apply(beta_hat, 1, \(th) t(th) %*% cov_X %*% th) - var_uncert
mean(V_I_hat)       
coda::HPDinterval(as.mcmc(V_I_hat)) 
var_uncert / mean(V_I_hat)


# We can do the same separately for the slope
V_S_hat <- (var(design_mat[ , "Nutrient"]))
pi  <- V_S_hat / V_I_hat
plot(as.mcmc(pi))
mean(pi)                 # 0.24
coda::HPDinterval(as.mcmc(pi)) # 0.18 - 0.30


beta_hat <- fixef(m1.brms, summary = FALSE)
# Setting up the posterior distribution as a list, it'll be convenient later
beta_hat <- apply(beta_hat, 1, \(vec) { vec }, simplify = FALSE)

# Getting the G-matrix
G_beta_hat <- VarCorr(m1.brms, summary = FALSE)[["Community"]][["cov"]] %>%
  apply(1, \(mat) { mat }, simplify = FALSE)




V_R <- m1.sum$sigma^2
V_E <- m1.sum$vcov[1,1]
V_S <- m1.sum$vcov[2,2]
cov_ES <- m1.sum$vcov[2,1]
R_ES <- cov_ES / sqrt(V_E * V_S) 
R_ES


V_E / (V_E + V_S + 2*cov_ES)
V_S / (V_E + V_S + 2*cov_ES)
