setwd("/home/sur/lab/exp/2026/today/")
library(tidyverse)
library(patchwork)
library(brms)
library(Reacnorm)
box::use(./functions/crn)


#' # Simulate some simple data
#' A bit cumbersome, but the line below simulate data for 4 **id**'s
#' (could be strains or syncoms) in 3 temps, 2 reps, and 4 batches (1 rep
#' of 2 id's per batch). A switch in id's per batch is simulated.
#' Overall variability is modelled (includes measurement error) and
#' batch-to-batch variations is simulated as additive 4-times smaller
#' than replicate variability (in reality it could be opposite, more batch
#' than biological variability)
set.seed(1234567)
var_rep <- 0.2
var_batch <- 0.05
b_batch <- rnorm(n = 4, mean = 0, sd = sqrt(var_batch))
names(b_batch) <- c("A","B","C","D")
Dat <- bind_rows(
  # First rep
  tibble(id = "ST01",
         temp = c(30,37,42),
         AUC = c(4,4,4),
         batch = "A"),
  tibble(id = "ST02",
         temp = c(30,37,42),
         AUC = c(4,4.5,5),
         batch = "A"),
  tibble(id = "ST03",
         temp = c(30,37,42),
         AUC = c(4,4.5,3.7),
         batch = "B"),
  
  tibble(id = "ST04",
         temp = c(30,37,42),
         AUC = c(4,4.1,3.2),
         batch = "B"),
  
  # Second rep with flipped batches
  tibble(id = "ST01",
         temp = c(30,37,42),
         AUC = c(4,4,4),
         batch = "C"),
  tibble(id = "ST02",
         temp = c(30,37,42),
         AUC = c(4,4.5,5),
         batch = "D"),
  tibble(id = "ST03",
         temp = c(30,37,42),
         AUC = c(4,4.5,3.7),
         batch = "C"),
  tibble(id = "ST04",
         temp = c(30,37,42),
         AUC = c(4,4.1,3.2),
         batch = "D")
) %>%
  mutate(AUC = AUC + rnorm(n = length(id), mean = 0, sd = var_rep)) %>%
  mutate(AUC = AUC + b_batch[batch])

#' Contigency table of samples
ftable(id + batch ~ temp, data = Dat)

#' Plot of the resulting *reaction norms*.
p1 <- Dat %>%
  ggplot(aes(x = temp, y = AUC,
             group = interaction(id, batch, sep = "_", drop = TRUE))) +
  geom_line(aes(col = id), linewidth = 3) +
  theme_classic()
p1

#' # Reaction norm
#' First it is convenient to convert temp to -1, 0, 1 for low, med and high
#' TO DO: test if normalizing AUC also helps

Dat <- Dat %>%
  mutate(temp = replace(temp, temp == 30, -1)) %>%
  mutate(temp = replace(temp, temp == 37, 0)) %>%
  mutate(temp = replace(temp, temp == 42, 1))
Dat

#' Try some general models:
#' 
#' * model 0: A simple mopdel where the effect of temperature is constant
#' for all id's with random intercepts for each id and batch
#' * model 1: A model where there is different polynomial effect of temperature
#' on each id, there is a separate batch random intercept.
#' * model 2: A model with an overall temperature effect, and a id specific
#' intercept and temperature effect. Temperature effects are polynomial of
#' order 1. There is a separate batch effect
#' * model 3: Same as model 2 but with polynomial effects of order 2. Allows
#' for curvature in reaction norm
#' 
#' Models 2 & 3 are models *a la Villemereuil* and can be easily partitioned
#' into different sources of variance. I expect model 3 will be the most
#' appropriate for our data. Note that in order to use model 3 we need
#' to have at least 3 environmental variables.
model_f0 <- formula(AUC ~ 1 + temp + temp^2 + (1 | id) + (1 | batch))
model_f1 <- formula(AUC ~ 1 + (1 + temp + temp^2 | id) + (1 | batch))
model_f2 <- formula(AUC ~ 1 + temp + (1 + temp | id) + (1 | batch))
model_f3 <- formula(AUC ~ 1 + temp + temp^2 + (1 + temp + temp^2 | id) + (1 | batch))


mp0 <- brm(model_f0,
           data = Dat,
           chains = 4, iter = 4000, warmup = 3000, cores = 4,
           control = list(adapt_delta = 0.99),
           save_pars = save_pars(all = TRUE))
mp1 <- brm(model_f1,
           data = Dat,
           chains = 4, iter = 4000, warmup = 3000, cores = 4,
           control = list(adapt_delta = 0.99),
           save_pars = save_pars(all = TRUE))
mp2 <- brm(model_f2,
           data = Dat,
           chains = 4, iter = 4000, warmup = 3000, cores = 4,
           control = list(adapt_delta = 0.99),
           save_pars = save_pars(all = TRUE))
mp3 <- brm(model_f3,
           data = Dat,
           chains = 4, iter = 4000, warmup = 3000, cores = 4,
           control = list(adapt_delta = 0.99),
           save_pars = save_pars(all = TRUE))

# summary(mp0)
# summary(mp1)
# summary(mp2)
# summary(mp3)

#' As expected mp3 is the best model, though it is tied with mp1 in the simulated
#' data
LOO(mp0, mp1, mp2, mp3, moment_match = TRUE, reloo = TRUE)

#' In any case simpler don't capture behavior as expected. I will
#' calculate the contributions of different factors for models 1 & 3. In
#' real data we only need to calculate whatever is the best model...unless
#' there is no clear best



# crn$partition_variance_polynomial(mp = mp1, pheno_name = "AUC", com_name = "id")

# Get formula for design i
design_f <- reformulas::findbars(model_f1)[[1]][[2]]
design_f


#' # Using reacnorm package

#' We need to define a matrix of relatedness between id's (strains or syncoms
#' depending on dataset). For strains we would use the gANI, for syncoms we
#' can use the proportion of shared strains or the UniFrac distance. 
#' Here, based on the values chose for simulation I will assume that ST04 & ST03
#' are quite similar, ST01 is a bit more disimilar to both, and ST01 is the
#' most different of all,
A <- diag(1, nrow = length(unique(Dat$id)))
colnames(A) <- rownames(A) <- unique(Dat$id)
A["ST04", "ST03"] <- A["ST03", "ST04"] <- 0.9
A["ST04", "ST01"] <- A["ST01", "ST04"] <- 0.7
A["ST03", "ST01"] <- A["ST01", "ST03"] <- 0.7
A["ST02", "ST01"] <- A["ST01", "ST02"] <- 0.5
A["ST02", "ST03"] <- A["ST03", "ST02"] <- 0.5
A["ST02", "ST04"] <- A["ST04", "ST02"] <- 0.5
A

#' Since we are going to use a quadratic model (to incorporate curvature),
#' it is convenient to square temperature
Dat <- Dat %>%
  mutate(temp_sq = temp ^ 2)

#' Define model formula. IMPORTANT: Incorporate our relatedness matrix
#' in the grouping factor. IMPORTANT:2  the slope temrs (tmp & temp_sq here),
#' should be identical in the fixed and random effects
model_f <- brmsformula(AUC ~ 1 + temp + temp_sq + ( 1 + temp + temp_sq | gr(id, cov = A) ))
m.quad_crn <- brm(model_f,
                  data = Dat,
                  data2 = list(A = A),
                  save_pars = save_pars(group = FALSE),
                  chains = 4,
                  cores = 4,
                  seed = 6543,
                  iter = 5000,
                  warmup = 3000,
                  control = list(adapt_delta = 0.99))
summary(m.quad_crn)
plot(m.quad_crn)

#' Same but with batch effects
model_f <- brmsformula(AUC ~ 1 + temp + temp_sq + ( 1 + temp + temp_sq | gr(id, cov = A) ) + (1 | batch))
m.quad_crn_batch <- brm(model_f,
                  data = Dat,
                  data2 = list(A = A),
                  save_pars = save_pars(group = FALSE),
                  chains = 4,
                  cores = 4,
                  seed = 6543,
                  iter = 5000,
                  warmup = 3000,
                  control = list(adapt_delta = 0.99))

summary(m.quad_crn_batch)
plot(m.quad_crn_batch)

#' In this case the warnings about treedepth probably have to do
#' with the fact that the relatedness matrix makes no sense with the observations
#' in real data we need to pay attention to warnings. I'll ignore here.

#' Compare models
# LOO(m.quad_crn, m.quad_crn_batch)

#' In real data we are likely to have batch effects so I will focus on the model
#' with batch term


#' # Plotting the reaction norm
#' First we calculate the 95% posterior intervals
Preds <- Dat %>%
  mutate(preds = predict(m.quad_crn_batch, re_formula = NA) %>%
           as_tibble()) %>%
  tidyr::unpack(preds) %>%
  select(temp,
         preds = Estimate,
         preds_low = Q2.5,
         preds_up = Q97.5) %>%
  summarise(across(starts_with("preds"), mean),
            .by = temp)
Preds

#' Then we add to our base plot
p1 <- Dat %>%
  ggplot(aes(x = temp, y = AUC)) +
  geom_line(aes(col = id,
                group = interaction(id, batch, sep = "_", drop = TRUE)),
            linewidth = 3) +
  geom_ribbon(data = Preds,
              mapping = aes(x = temp, ymin = preds_low, ymax = preds_up, y = preds),
              alpha = 0.3) +
  geom_line(data = Preds,
            mapping = aes(x = temp, y = preds),
            linewidth = 1) + 
  
  theme_classic()
p1

#' Decompose the variance using the Reacnorm package
seq_env <- c(-1, 0,1)
# seq_env <- seq(from = -1, to = 1, by = 0.1)
env_X <- cbind(1, seq_env, seq_env ^ 2) # Design matrix for the quadratic model

#' We extract the different models parameters. IMPORTANT: We need to change
#' the names of the parameters to "a" (for intercept), "b" (for the main env
#' variable), and "c" (for quadratic term of the environment). The order will
#' match the order of terms in the formula of the model fit

#' Firs the fixed effect parameters
theta <- fixef(m.quad_crn_batch, robust = TRUE)[, "Estimate"] # Median estimates
names(theta) <- c("a", "b", "c")
theta               

#' We get the estimate uncertainties as well.
theta_vcov <- vcov(m.quad_crn_batch)
rownames(theta_vcov) <- colnames(theta_vcov) <- names(theta)
theta_vcov

#' The matrix of relatedness (G-matrix) of the parameters
G_mat <-  VarCorr(m.quad_crn_batch, robust = TRUE)[["id"]][["cov"]][ , "Estimate", ]
rownames(G_mat) <- colnames(theta_vcov) <- names(theta)
G_mat

#' The residual and batch SD's which we square (for variance) and add (for total variance)
vr_ext <- VarCorr(m.quad_crn_batch, robust = TRUE)[["residual__"]][["sd"]][ , "Estimate" ] ^2 +
  VarCorr(m.quad_crn_batch, robust = TRUE)[["batch"]][["sd"]][,"Estimate"] ^ 2
vr_ext

#' Decompose the variance. The `wt_env` parameter is designed for natural
#' env distributions. Here, since it is an experiment, and we don't know
#' the natural distributions, we give equal weights to all envs
vplas <- rn_pi_decomp(theta = theta,
                      V_theta = G_mat,
                      env = seq_env,
                      shape = expression(a + b * x + c * x^2),
                      wt_env = rep(1, times = length(seq_env)))
vplas

#' Here we see that there is little overall variation due to the environment (plasticity),
#' around 1.5% (V_plas), this makes sense as the variation between temps, is much smaller
#' than the variation between ids, in real data it could be quite different.
#' Then, Pi_Sl is the proportion of V_Plas is explained, and Pi_Cv is the
#' proportion of V_Plas explained by the curvature. Here There is no common slope
#' so that is why almost 70% of V_Plas is explained by the curvature.
#' 
#' Thecnically we could also use the  Phi decomposition to reach a similar
#' conclusion. ¿Or only if wt_env is normal?
# rn_phi_decomp(theta = theta,
#               X = env_X,
#               S = theta_vcov,
#               wt_env = rep(1, times = length(seq_env)))

#'  # Relatedness decomposition
vrel <- rn_gen_decomp(theta = theta,
                      G_theta = G_mat,
                      X = seq_X,
                      wt_env = rep(1, times = length(seq_env)))
vrel

#' Here V_Add is the the variance due to differences between ids (here strains).
#' Which can be decomposed as V_A, the variance due to difference between mean
#' phenotypic values of each id (also called environment-blind), and V_AxE
#' which is the variance around those means (the difference in the plastic response
#' to environment between ids). Can be expressed as a percentace, and here
#' we would see that about a third (~32%) of the variance between id's is due
#' to differences in their response to the environment, and the remaining, is
#' do to overall differences in their mean phenotypic values.
#' 
#' The gamma and iota values further decompose V_A & V_AxE, respectively, 
#' into their slope (Gamma_b, Iota_b) and curvature (Gamma_c, Iota_c) components.
#' and curvature elements. Though negative values have to be treated with care.

#' We can normalize everything as a function of the total phenotypic variance
#' (including batch and residual)
var_tot <- vplas[["V_Plas"]] + vrel[["V_Add"]] + vr_ext
var_pheno <-
  c(P2 = vplas[["V_Plas"]] / var_tot,
    h2_RN = vrel[["V_Add"]] / var_tot,
    h2 = vrel[["V_A"]] / var_tot,
    h2_I = vrel[["V_AxE"]] / var_tot,
    T2 = (vplas[["V_Plas"]] + vrel[["V_Add"]]) / var_tot)
var_pheno

#' Getting the posterior

#' Firs the fixed effect parameters
theta_post <- fixef(m.quad_crn_batch, summary = FALSE)
colnames(theta_post) <- c("a", "b", "c")
head(theta_post)

#' The residual and batch SD's (squared and added)
vr_ext_post <- VarCorr(m.quad_crn_batch, summary = FALSE)[["residual__"]][["sd"]][ , 1 ] ^ 2 +
  VarCorr(m.quad_crn_batch, summary = FALSE)[["batch"]][["sd"]][ , 1 ] ^ 2
head(vr_ext_post)

#' G_mat needs to be transformed into a list
G_mat_post <- VarCorr(m.quad_crn_batch, summary = FALSE)[["id"]][["cov"]] %>%
  apply(1,function(mat){mat}, simplify = FALSE) %>% # Converts 3D array into list
  map(function(mat){
    rownames(mat) <- colnames(mat) <- c("a", "b", "c")
    return(mat)
    })
head(G_mat_post)

#' For convenience, combine everything into a posterior distribution object
#' using the posterior package
Post <- as_draws_df(theta_post)
Post[["G"]] <- G_mat_post
Post[["V_R"]] <- vr_ext_post
Post <- posterior::thin_draws(Post, thin = nrow(theta_post) / 1000)
post_info <- select(Post, starts_with(".")) # convenience for new objects
Post


#' VPlas decomposition
vplas_post <- Post %>%
  pmap(function(a, b, c, G, V_R, .chain, .iteration, .draw){
    # rn_phi_decomp(theta = c(a = a, b = b, c = c),
    #               X = seq_X,
    #               S = theta_vcov,
    #               wt_env = rep(1, times = length(seq_env)))
    
    rn_pi_decomp(theta = c(a = a, b = b, c = c),
                 V_theta = G,
                 env = seq_env,
                 shape = expression(a + b * x + c * x^2),
                 wt_env = rep(1, times = length(seq_env)))
    }, .progress = TRUE) %>%
  bind_rows() %>%
  select(where(function(column){abs(mean(column)) > 1e-5})) %>%
  cbind(post_info) %>% # Add chain and draw info
  as_draws_df()
vplas_post


posterior::summarise_draws(vplas_post)
bayesplot::mcmc_trace(vplas_post)

bayesplot::mcmc_areas(vplas_post,
                      regex_pars = "^V",
                      prob = 0.95,
                      area_method = "scaled height") /
bayesplot::mcmc_areas(vplas_post,
             regex_pars = "^[^V]",
             prob = 0.95,
             area_method = "scaled height") +
  patchwork::plot_layout(heights = c(1, 2))





