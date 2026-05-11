setwd("~/lab/exp/2026/today/")
library(tidyverse)
library(miaSim)

set.seed(654)
n_sim <- 10
n_spec <- 10
n_passages <- 3
passage_dil <- 1e3
t_passage <- 100
Syncoms <- expand_grid(S1 = c(1,6),
                       S2 = c(2,7),
                       S3 = c(3,8),
                       S4 = c(4,9),
                       S5 = c(5,10))
Syncoms
M <- randomA(n_species = n_spec,
             diagonal = -0.5,
             connectance = 0.5)

date()
Sims <- NULL
for(sc in 1:nrow(Syncoms)){
  # sc <- 1
  cat("=== SynCom ", sc, "\n")
  specs <- Syncoms[sc,] %>% unlist
  M_sc <- M[specs, specs]
  abuns_p <- rep(0.001, nrow(M_sc))
  hrs <- 0
  res <- tibble(time = hrs, spec = row.names(M_sc), abun = abuns_p)
  for(p in 1:n_passages){
    Sim <- simulateGLV(n_species = nrow(M_sc),
                       names_species = row.names(M_sc),
                       A = M_sc,
                       x0 = abuns_p,
                       t_end = t_passage)
    abuns_p <- assay(Sim)[,1000]
    hrs <- hrs + t_passage
    res <- bind_rows(res,
                     tibble(time = hrs, 
                            spec = row.names(M_sc),
                            abun = abuns_p))
    abuns_p <- abuns_p / passage_dil
  }
  res$syncom <- sc
  Sims <- bind_rows(Sims, res)
  res <- NULL
}
Sims
date()



strain_list <- 1:10
selected_time <- 300


Res <- NULL 
for(s1 in strain_list){
  for(s2 in strain_list){
    # s1 <- 1
    # s2 <- 2
    if(s1 == s2){
      next
    }
    cat(s1, s2, "\n")
    s1_ii <- apply(Syncoms, 1, function(x){s1 %in% x})
    s2_ii <- apply(Syncoms, 1, function(x){s2 %in% x})
    
    s1_sc <- which(s1_ii)
    s2_sc <- which(s2_ii)
    # s1_sc <- colnames(Syncoms)[ which(Syncoms[s1_ii,] == 1) ]
    # s2_sc <- colnames(Syncoms)[ which(Syncoms[s2_ii,] == 1) ]
    
    int_sc <- intersect(s1_sc, s2_sc)
    s1_only <- setdiff(s1_sc, int_sc)
    
    if(length(int_sc) && length(s1_only)){
      # Sam <- Meta %>%
      #   filter(hrs %in% selected_time) %>%
      #   filter(Community %in% union(s1_only, int_sc))
      
      Sam <- Sims %>%
        filter(syncom %in% union(s1_only, int_sc)) %>%
        filter(time %in% selected_time) %>%
        mutate(pair = 1*(syncom %in% int_sc))
      Sam
      
      # Test <- Freqs %>%
      #   filter(Strain %in% c(s1)) %>%
      #   filter(label %in% Sam$label) %>%
      #   left_join(Sam %>%
      #               select(label, sample, Community, temp, hrs, exp), by = "label") %>%
      #   mutate(pair = 1*(Community %in% int_sc))
      # Test %>%
      #   print(n = 1000)
      
      
      m1 <- lm(log2(abun + 1e-8) ~ pair + syncom , data = Sam )
      summary(m1)
      
      res <- tibble(s1 = s1,
                    s2 = s2,
                    log2FC = coef(m1)["pair"],
                    pval = summary(m1)$coefficients["pair",4],
                    int = M[s1, s2])
      
      # ggplot(Test,aes(x = factor(pair), y = abs1)) +
      #   geom_boxplot() +
      #   geom_point() +
      #   theme_classic()
      
      Res <- bind_rows(Res, res)
      
    }
    
    # res <- tibble(s1 = s1,
    #               s2 = s2,
    #               log2FC = NA,
    #               pval = NA,
    #               int = M[s1, s2])
    # 
    # Res <- bind_rows(Res, res)
    
  }
  
}
Res %>% print(n = 1000)
Res %>%
  ggplot(aes(x = log2FC, y = int)) +
  geom_point()
cor(Res$log2FC, Res$int)
