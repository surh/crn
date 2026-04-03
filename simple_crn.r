setwd("/home/sur/lab/exp/2026/today/")
library(tidyverse)
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
#' can use the proporion of shared strains or the UniFrac distance. Here for
#' simplicity I would assume that ids are all equally dissimilar (technically
#' I'm assuming iid). NOTE: Need to incorporate this with real data.
A <- diag(1, nrow = length(unique(Dat$id)))
colnames(A) <- rownames(A) <- unique(Dat$id)
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

#' Plotting the reaction norm

Preds <- Dat %>%
  mutate(preds = predict(m.quad_crn, re_formula = NA) %>%
           as_tibble()) %>%
  tidyr::unpack(preds) %>%
  select(temp,
         preds = Estimate,
         preds_low = Q2.5,
         preds_up = Q97.5) %>%
  summarise(across(starts_with("preds"), mean),
            .by = temp)
Preds

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


