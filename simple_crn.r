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
summary(m.quad_crn_batch)

#+ Traceplots
plot(m.quad_crn_batch)

#+ Select model
# Convenient in case there is more than one model
main_model <- m.quad_crn_batch
save(main_model, file = "main_model.rdat") # To avoid refitting

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
#' We need to extract some values from the model fit. For full bayesian
#' treatment we need the posterior estimates of each parameter. We have to add
#' the residual and batch variance, so we extract the SDs, square them, and
#' add them
#+ Extract model params
# Constructed values
seq_env <- c(-1, 0,1)
env_X <- cbind(1, seq_env, seq_env ^ 2) # Design matrix for the quadratic model
theta_vcov <- vcov(main_model) # Uncertainty in fixed params

# Fixed effect posterior
theta_post <- fixef(main_model, summary = FALSE)
colnames(theta_post) <- c("a", "b", "c") # Rename for phi decomposition
head(theta_post)

# G_mat needs to be transformed from a 3D array into a list
G_mat_post <- VarCorr(main_model, summary = FALSE)[["id"]][["cov"]] %>%
  apply(1,function(mat){mat}, simplify = FALSE) %>% # Converts 3D array into list
  map(function(mat){
    rownames(mat) <- colnames(mat) <- c("a", "b", "c") # Rename for phi-decomp
    return(mat)
  })
head(G_mat_post)

# Batch and residual variance (aka external variance)
var_ext_post <- VarCorr(main_model, summary = FALSE)[["residual__"]][["sd"]][ , 1 ] ^ 2 +
  VarCorr(main_model, summary = FALSE)[["batch"]][["sd"]][ , 1 ] ^ 2
head(var_ext_post)

#' Now we combine everything into one big posterior object. This requires
#' the `posterior` package. It simplifies keeping chains and iterations
#' for plotting with `bayesplot`.
#+ Overall posterior
Post <- as_draws_df(theta_post)
Post[["G"]] <- G_mat_post
Post[["V_R"]] <- var_ext_post
# Post <- posterior::thin_draws(Post, thin = nrow(theta_post) / 250) # For debugging
post_info <- select(Post, starts_with(".")) # convenience for new objects
Post

#' We decompose the mean reaction norm (Vplas). We could use the pi-decomposition
#' but the phi-decomposition is more general, and incorporates uncertainty of estimates
#+ RN decomposition on posterior
Vplas_post <- Post %>%
  pmap(function(a, b, c, G, V_R, .chain, .iteration, .draw, env_X, theta_vcov){
    rn_phi_decomp(theta = c(a = a, b = b, c = c),
                  X = env_X,
                  S = theta_vcov,
                  wt_env = rep(1, times = nrow(env_X)))
    
    # Pi decomposition is much slower, but seems to work better
    # in simulated data
    # rn_pi_decomp(theta = c(a = a, b = b, c = c),
    #              V_theta = G,
    #              env = seq_env,
    #              shape = expression(a + b * x + c * x^2),
    #              wt_env = rep(1, times = length(seq_env)))
  }, theta_vcov = theta_vcov, env_X = env_X, .progress = TRUE) %>%
  bind_rows() %>%
  cbind(post_info) %>% # Add chain and draw info
  as_draws_df() # Convert to posterior
Vplas_post

#' We repeat the iteration but now we do the relatedness ("gen") decomposition
#' to decompose the variation in plasticity between ids
#+ Variation in RN decomposition
Vrel_post <- Post %>%
  pmap(function(a, b, c, G, V_R, .chain, .iteration, .draw, env_X){
    rn_gen_decomp(theta = c(a = a, b = b, c = c),
                  G_theta = G,
                  X = env_X,
                  wt_env = rep(1, times = nrow(env_X)))
  }, env_X = env_X, .progress = TRUE) %>%
  bind_rows() %>%
  cbind(post_info) %>% # Add chain and draw info
  as_draws_df() # Convert to posterior
Vrel_post

#' Now we look at the results, for `Vplas_post` we are mostly interested
#' in V_Plas, Phi_b and Phi_c, and we expect Phi_b_c to be very small. We can
#' check that with posterior::sumarise and plot only the desired variables.
#' Because V_Plas is by definition in a different scale, we plot it separatedly.
#+ Check Vplas
posterior::summarise_draws(Vplas_post)
bayesplot::mcmc_trace(Vplas_post,
                      pars = c("V_Plas", "Phi_b", "Phi_c"))
bayesplot::mcmc_areas(Vplas_post,
                      pars = c("V_Plas", "Phi_b", "Phi_c"),
                      prob = 0.8,
                      prob_outer = 0.9,
                      point_est = "median",
                      area_method = "equal area")

#' This fit is problematic because it gives us a negative VPlas (pi-decomp doesn't),
#' but the imterpretation of Phi_b and Phi_c is that variation in the average
#' reaction norm is greater due to the slope (Phi_b) than to the curvature (Phi_c).
#' That said, the model has problems here because the average reaction norm is almost
#' a straight horizontal line

#' Then we check the variation between ids. There are potentially many
#' parameters here. So here we select the variables that have an absolute
#' median posterior value above 1e-3 to docus on that. Also because V_Add,
#' V_A and V_AxE can often be in a different range than the rest of the parameters
#' we plot them by themselves
#+ Check Vrel
posterior::summarise_draws(Vrel_post)
vars_to_plot <- Vrel_post %>%
  select(where(function(column){abs(median(column)) > 1e-3})) %>%
  select(!starts_with(".")) %>% colnames()

bayesplot::mcmc_trace(Vrel_post,
                      pars = vars_to_plot)
bayesplot::mcmc_areas(Vrel_post,
                      pars = c("V_Add", "V_A", "V_AxE"),
                      prob = 0.8,
                      prob_outer = 0.9,
                      point_est = "median",
                      area_method = "equal area")
bayesplot::mcmc_areas(Vrel_post,
                      pars = setdiff(vars_to_plot, c("V_Add", "V_A", "V_AxE")),
                      prob = 0.8,
                      prob_outer = 0.9,
                      point_est = "median",
                      area_method = "equal area")
#' Basically we see that the there is a significant amount of variation
#' explained by differences between strains (V_Add), and that more of that
#' variation can be explained by difference between the average phenotypes of
#' each strain (V_A), and a smaller fraction because of differences in the change
#' of the phenotypes (response) with respect to the environment between strains
#' (V_AxE). Then each of those can be further decomposed (gammas and iotas).

#' Normally we want to express the main variance components as proportion of
#' the total variance, so wee need to add the different types of variance and
#' normalize everything. We use the posterior package. Some of the standardized
#' values have standard names in popgen theory (e.g. V_A / V_Tot = H^2 = heritability)
#+ Standardize variance
Var_std_post <- posterior::bind_draws(Post, Vplas_post, Vrel_post) %>%
  posterior::subset_draws(c("V_Plas", "V_Add", "V_A", "V_AxE", "V_R")) %>%
  posterior::mutate_variables(V_Tot = V_Plas + V_Add + V_R) %>%
  transmute(P2 = V_Plas / V_Tot,
            H2_RN = V_Add / V_Tot,
            H2 = V_A / V_Tot,
            H2_I = V_AxE / V_Tot,
            T2 = (V_Plas + V_Add) / V_Tot) %>%
  cbind(post_info) %>% # Add draws info
  as_draws_df() # Cionvert to posterior

#' We summarise and plot the results.
#+ Inspect standardized values
posterior::summarise_draws(Var_std_post)
bayesplot::mcmc_trace(Var_std_post)
bayesplot::mcmc_areas(Var_std_post,
                      prob = 0.8,
                      prob_outer = 0.9,
                      point_est = "median",
                      area_method = "equal area")

#' Here we get very bad values, with numbers above 1 and below 0 (they are
#' supposed to be proportions), thisi is because V_Plas is negative. In real
#' data it shouldn't be negative.
date()
