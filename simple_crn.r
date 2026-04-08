#+ Load dependencies
library(tidyverse)
library(brms)
library(Reacnorm)
date()

#' # Simulate some simple data
#' A bit cumbersome, but the line below simulate data for 4 **id**'s
#' (could be strains or syncoms) in 3 temps, 2 reps, and 4 batches (1 rep
#' of 2 id's per batch). A switch in id's per batch is simulated.
#' Overall variability is modeled (includes measurement error) and
#' batch-to-batch variations is simulated as additive 4-times smaller
#' than replicate variability (in reality it could be opposite, more batch
#' than biological variability)
#+ Simulate data
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

#' Plot of the resulting reaction norms
#+ Plot data
p1 <- Dat %>%
  ggplot(aes(x = temp, y = AUC,
             group = interaction(id, batch, sep = "_", drop = TRUE))) +
  geom_line(aes(col = id), linewidth = 3) +
  theme_classic()
p1

#' # Reaction norms
#' 
#' ## Prepare data
#' First it is convenient to convert temp to -1, 0, 1 for low, med and high.
#' 
#' **TO DO**: We should try the same analysis but standardizing the AUC values
#' as well.
#+ reformat data
Dat <- Dat %>%
  mutate(temp = replace(temp, temp == 30, -1)) %>%
  mutate(temp = replace(temp, temp == 37, 0)) %>%
  mutate(temp = replace(temp, temp == 42, 1))
Dat

#' ## Prepare variables
#' 
#' We need to define a matrix of relatedness between id's (strains or syncoms
#' depending on dataset). For strains we should try the following:
#' 1. The gANI
#' 2. Proportion of shared genes
#' 
#' For syncoms we have several options to test:
#' 1. The proportion of shared strains.
#' 2. The weighted UniFrac distance.
#' 3. The Bray-Curtis dissimilarity between inocula
#' 
#' Here, I make up some nombers. I will assume that ST04 & ST03
#' are quite similar, ST01 is a bit more dissimilar to both, and ST01 is the
#' most different of all. **IMPORTANT**: Because the numbers are made up,
#' the model don't fit very well and I get a number of warnings, but we
#' expect inr eal data we would get a better fit.
#' 
#' **QUESTION** Does the matrix need to be symetrical?
#+ Relatedness matrix
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
#' it is convenient to have a separate variable with square temperature-
#+ Square temperature
Dat <- Dat %>%
  mutate(temp_sq = temp ^ 2)

#' ## Run models
#' 
#' Just one model is tried here, a quadratic model (temp + temp^2), with
#' a random intercept for batch effects. There is a bunch of warnings here
#' and in the subsequent analysis. I ignore them because the data is simulated
#' and so I don't expect this to be the right model, but with real data we
#' have to take the warnings seriously
#+ Fit model
model_f <- brmsformula(AUC ~ 1 + temp + temp_sq + ( 1 + temp + temp_sq | gr(id, cov = A) ) + (1 | batch))
m.quad_crn_batch <- brm(model_f,
                  data = Dat,
                  data2 = list(A = A),
                  save_pars = save_pars(all = TRUE),
                  chains = 4,
                  cores = 4,
                  seed = 6543,
                  iter = 5000,
                  warmup = 3000,
                  control = list(adapt_delta = 0.99))

#+ Model summaries
summary(m.quad_crn)

#+ Traceplots
plot(m.quad_crn)
plot(m.quad_crn_batch)

#+ Select model
# Convenient in case there is more than one model
main_model <- m.quad_crn_batch

#' ## Plot average reaction norm
#' First we calculate the 95% posterior intervals and then we add them to
#' our base plot
#+ Plot average RN
Preds <- Dat %>%
  mutate(preds = predict(main_model, re_formula = NA) %>%
           as_tibble()) %>%
  tidyr::unpack(preds) %>%
  select(temp,
         preds = Estimate,
         preds_low = Q2.5,
         preds_up = Q97.5) %>%
  summarise(across(starts_with("preds"), mean),
            .by = temp)

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

#' ## Decompose the variance using the Reacnorm package
#' We need to extract some values from the model fit
#+ Extract model params
#+ Env values
seq_env <- c(-1, 0,1)
env_X <- cbind(1, seq_env, seq_env ^ 2) # Design matrix for the quadratic model

#+ Fixed effect estimates Extract central estimates
theta <- fixef(main_model, robust = TRUE)[, "Estimate"] # Median estimates
names(theta) <- c("a", "b", "c") #' Only for polynomial degree 2
theta_vcov <- vcov(main_model)
rownames(theta_vcov) <- colnames(theta_vcov) <- names(theta)

#+ G-matrix
G_mat <-  VarCorr(main_model, robust = TRUE)[["id"]][["cov"]][ , "Estimate", ]
rownames(G_mat) <- colnames(G_mat) <- names(theta) # For polynomial random effects

#+ Residual and batch SD's
vr_ext <- VarCorr(main_model, robust = TRUE)[["residual__"]][["sd"]][ , "Estimate" ] ^ 2 +
  VarCorr(main_model, robust = TRUE)[["batch"]][["sd"]][,"Estimate"] ^ 2

#' Decompose the variance. The `wt_env` parameter is designed for natural
#' env distributions. Here, since it is an experiment, and we don't know
#' the natural distributions, we give equal weights to all envs
#+ Decompose variance
vplas <- rn_pi_decomp(theta = theta,
                      V_theta = G_mat,
                      env = seq_env,
                      shape = expression(a + b * x + c * x^2),
                      # shape = expression(a), # for polynomial 0
                      wt_env = rep(1, times = length(seq_env)))
vplas

# m.quad_crn pi
# V_Plas     Pi_Sl     Pi_Cv
# 1 0.02047858 0.2078569 0.7894149

# m.quad_crn_batch pi
# V_Plas     Pi_Sl     Pi_Cv
# 1 0.01522686 0.3046904 0.6946678

# m.p0_crn_batch
# V_Plas Pi_Sl Pi_Cv
# 1      0   NaN   NaN

#' Here we see that there is little overall variation due to the environment (plasticity),
#' around 1.5% (V_plas), this makes sense looking at the plot. Then, Pi_Sl
#' is the proportion of V_Plas is explained, and Pi_Cv is the
#' proportion of V_Plas explained by the curvature. Here There is no common slope
#' so that is why almost 70% of V_Plas is explained by the curvature.
#' 
#' Thecnically we could also use the  Phi decomposition to reach a similar
#' conclusion. ¿Or only if wt_env is normal?
rn_phi_decomp(theta = theta,
              X = env_X,
              S = theta_vcov,
              wt_env = rep(1, times = length(seq_env)))

# m.quad_crn
# V_Plas     Phi_b     Phi_c Phi_b_c
# 1 -0.5871767 0.7613559 0.2386441       0
# 
# m.quad_crn_batch
# V_Plas     Phi_b     Phi_c Phi_b_c
# 1 -0.6560756 0.7010698 0.2989302       0



#'  # Relatedness decomposition
#+ Relatedness decomposition
vrel <- rn_gen_decomp(theta = theta,
                      G_theta = G_mat,
                      X = env_X,
                      wt_env = rep(1, times = length(seq_env)))
vrel

# m.quad_crn 
# V_Add       V_A     V_AxE   Gamma_a   Gamma_b   Gamma_c Gamma_a_b  Gamma_a_c Gamma_b_c Iota_a
# 1 1.290644 0.8594782 0.4311653 0.5370563 0.2398479 0.2826664         0 -0.0595706         0      0
# Iota_b   Iota_c Iota_a_b Iota_a_c Iota_b_c
# 1 0.717957 0.282043        0        0        0

# m.quad_crn_batch
# V_Add      V_A     V_AxE   Gamma_a   Gamma_b   Gamma_c Gamma_a_b  Gamma_a_c Gamma_b_c Iota_a
# 1 1.628748 1.112221 0.5165265 0.6347726 0.2013925 0.3472158         0 -0.1833808         0      0
# Iota_b   Iota_c Iota_a_b Iota_a_c Iota_b_c
# 1 0.635045 0.364955        0        0        0

m.p0_crn_batch

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
#+ Standardize variance
var_tot <- vplas[["V_Plas"]] + vrel[["V_Add"]] + vr_ext
var_pheno <-
  c(P2 = vplas[["V_Plas"]] / var_tot,
    h2_RN = vrel[["V_Add"]] / var_tot,
    h2 = vrel[["V_A"]] / var_tot,
    h2_I = vrel[["V_AxE"]] / var_tot,
    T2 = (vplas[["V_Plas"]] + vrel[["V_Add"]]) / var_tot)
var_pheno

#' Getting the posterior

#' First the fixed effect parameters
#+ Get posterior distributions
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
#+ RN decomposition on posterior
vplas_post <- Post %>%
  pmap(function(a, b, c, G, V_R, .chain, .iteration, .draw){
    rn_phi_decomp(theta = c(a = a, b = b, c = c),
                  X = env_X,
                  S = theta_vcov,
                  wt_env = rep(1, times = length(seq_env)))

    # Pi decomposition is much slower, but seems to work better
    # in simulated data
    # rn_pi_decomp(theta = c(a = a, b = b, c = c),
    #              V_theta = G,
    #              env = seq_env,
    #              shape = expression(a + b * x + c * x^2),
    #              wt_env = rep(1, times = length(seq_env)))
    }, .progress = TRUE) %>%
  bind_rows() %>%
  select(where(function(column){abs(mean(column)) > 1e-5})) %>%
  cbind(post_info) %>% # Add chain and draw info
  as_draws_df()
vplas_post

#+ Plot posterior distribution of decomp
posterior::summarise_draws(vplas_post)
bayesplot::mcmc_trace(vplas_post)

bayesplot::mcmc_areas(vplas_post,
                      pars = "V_Plas",
                      prob = 0.95,
                      area_method = "scaled height") /
bayesplot::mcmc_areas(vplas_post,
             pars = c("Phi_b", "Phi_c"),
             # pars = c("Pi_Sl", "Pi_Cv"),
             prob = 0.95,
             area_method = "scaled height") +
  patchwork::plot_layout(heights = c(1, 2))



#' Relatednes decomposition
#+ Posterior distribution of relatedness
vrel_post <- Post %>%
  pmap(function(a, b, c, G, V_R, .chain, .iteration, .draw){
    rn_gen_decomp(theta = c(a = a, b = b, c = c),
                  G_theta = G,
                  X = env_X,
                  wt_env = rep(1, times = length(seq_env)))
  }, .progress = TRUE) %>%
  bind_rows() %>%
  select(where(function(column){abs(mean(column)) > 1e-5})) %>%
  cbind(post_info) %>%
  as_draws_df()
vrel_post

posterior::summarise_draws(vrel_post)

#+ Plot of posterior distribution of relatedness
bayesplot::mcmc_trace(vrel_post)
bayesplot::mcmc_areas(vrel_post,
                      pars = c("V_Add", "V_A", "V_AxE"),
                      prob = 0.95,
                      area_method = "scaled height")
bayesplot::mcmc_areas(vrel_post,
                      regex_pars = "^[^V]",
                      prob = 0.95,
                      area_method = "scaled height") 

date()
