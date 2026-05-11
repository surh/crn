library(tidyverse)
setwd("~/lab/exp/2026/today/")

Freqs <- read_tsv("data/f_clean.tsv")
Freqs <- Freqs %>%
  rename(Strain = row.names) %>%
  mutate(Strain = replace(Strain, Strain == "NS_042g_27F", "ST00042")) %>%
  mutate(Strain = replace(Strain, Strain == "NS_164C_27F", "ST00164")) %>%
  mutate(Strain = replace(Strain, Strain == "NS_110C_1_27F", "ST00110"))
Freqs
  
Meta <- read_tsv("data/metadata_clean.tsv")
Meta <- Meta %>%
  mutate(sample = str_remove(label, "[AB]$")) %>%
  rename(Community = community)
Meta

Growth <- read_tsv("pilot_syncom_growth_curves.tsv")
Growth
Growth <- Growth %>%
  select(Community, temp, OD600, hrs = total_time_h, exp = batch)
Growth

Syncoms <- bind_cols(read_tsv("/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS1/syncoms.tsv") %>%
                       rename(Strain = strain),
                 read_tsv("/home/sur/lab/data/2024_rhizo_pilot_syncom_NS/NS2/syncoms.tsv")) 
Syncoms
Syncoms  <- Syncoms %>%
  select(-strain)
Syncoms



Freqs <- Freqs %>%
  pivot_longer(-Strain,
               names_to = "label",
               values_to = "count")
Meta <- Meta %>%
  left_join(Freqs %>%
              group_by(label) %>%
              summarise(depth = sum(count),
                        .groups = "drop"),
            by = "label")

Freqs <- Freqs %>%
  right_join(Meta %>% select(label, depth), by = "label") %>%
  mutate(freq = count / depth) %>%
  select(-depth)

Meta <- Meta %>%
  select(-day, -techrep) 
Meta


Growth <- Growth %>%
  filter(!(hrs == 0 & temp == 32)) %>%
  mutate(temp = replace(temp, hrs == 0, NA )) 
Growth %>% filter(hrs == 0) %>%
   print(n = 100)

Meta <- Meta %>%
  left_join(Growth, by = c("Community", "temp", "hrs", "exp"))

Freqs <- Freqs %>%
  left_join(Meta %>% select(label, OD600)) %>%
  mutate(abs1 = freq * OD600) %>%
  select(-OD600)




Syncoms[ is.na(Syncoms) ] <- 0
Syncoms

selected_time <- 48
strain_list <- Syncoms$Strain

Res <- NULL 
for(s1 in strain_list){
  for(s2 in strain_list){
    if(s1 == s2){
      next
    }
    cat(s1, s2, "\n")
    s1_ii <- which(Syncoms$Strain == s1)
    s2_ii <- which(Syncoms$Strain == s2)
    
    
    s1_sc <- colnames(Syncoms)[ which(Syncoms[s1_ii,] == 1) ]
    s2_sc <- colnames(Syncoms)[ which(Syncoms[s2_ii,] == 1) ]
    
    int_sc <- intersect(s1_sc, s2_sc)
    s1_only <- setdiff(s1_sc, int_sc)
    
    if(length(int_sc) && length(s1_only)){
      Sam <- Meta %>%
        filter(hrs %in% selected_time) %>%
        filter(Community %in% union(s1_only, int_sc))
      
      
      Test <- Freqs %>%
        filter(Strain %in% c(s1)) %>%
        filter(label %in% Sam$label) %>%
        left_join(Sam %>%
                    select(label, sample, Community, temp, hrs, exp), by = "label") %>%
        mutate(pair = 1*(Community %in% int_sc))
      # Test %>%
      #   print(n = 1000)
      
      
      m1 <- lm(log2(abs1 + 1e-8) ~ pair + Community + exp, data = Test )
      summary(m1)
      
      res <- tibble(s1 = s1,
                    s2 = s2,
                    log2FC = coef(m1)["pair"],
                    pval = summary(m1)$coefficients["pair",4])
      
      # ggplot(Test,aes(x = factor(pair), y = abs1)) +
      #   geom_boxplot() +
      #   geom_point() +
      #   theme_classic()
      
      Res <- bind_rows(Res, res)
      
    }
    
    res <- tibble(s1 = s1,
                  s2 = s2,
                  log2FC = NA,
                  pval = NA)
    
    Res <- bind_rows(Res, res)
    
  }
  
}
Res <- Res %>%
  filter(!is.na(log2FC)) %>%
  arrange(pval) %>%
  mutate(qval = p.adjust(pval, 'fdr')) %>%
  print(n = 100)
write_tsv(Res, "interactions.tsv")


Ints <- Res
p1 <- Ints %>%
  mutate(log2FC = replace(log2FC, qval > 0.05, 0)) %>%
  ggplot(aes(x = s1, y = s2)) +
  geom_tile(aes(fill = log2FC)) +
  scale_fill_gradient2() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank())
p1
ggsave("interaction_matrix.png", p1, width = 5, height = 4)


# mat <- Ints %>%
#   mutate(log2FC = replace(log2FC, qval > 0.05, 0)) %>%
#   select(-qval) %>%
#   # print(n = 1000) %>%
#   select(-pval) %>%
#   pivot_wider(names_from = s2, values_from = log2FC)

# row_strains <- mat$s1
# mat$s1 <- NULL
# mat <- as.matrix(mat)
# row.names(mat) <- row_strains
# mat <- t(mat)
# mat[is.na(mat)] <- 0
# mat
# 
# network <- igraph::graph_from_adjacency_matrix(mat, mode = "directed", weighted = TRUE)
# network
# plot(network)

library(igraph)
Mat <- Ints
Mat <- Mat %>%
  mutate(log2FC = replace(log2FC, qval > 0.05, 0)) %>%
  filter(abs(log2FC) > 0) %>%
  select(source = s2, target = s1, log2FC)
Mat
network <- graph_from_data_frame(Mat, directed = T)
co <- layout_with_fr(network,
                     minx = rep(-Inf, vcount(network)),
                     miny = rep(-Inf, vcount(network)),
                     maxx = rep(Inf, vcount(network)),
                     maxy = rep(Inf, vcount(network)))
png("network_all.png", width = 1000, height = 1000)
plot(network, layout = co, 
     edge.arrow.size = 2,
     edge.width = 4, edge.color = "black")
dev.off()

Mat %>% filter(target == "ST00046")

Nodestats <- NULL
Graphstats <- NULL
for(sc in colnames(Syncoms)[-1]){
  strains <- Syncoms$Strain[ which(Syncoms[[sc]] == 1) ]
  strains_ii <- which(V(network)$name %in% strains)
  n_com <- induced_subgraph(network, strains_ii)
  filename <- paste0("com_annot_", sc,"_network.png")
  png(filename, width = 1000, height = 1000)
  plot(n_com, layout = co[strains_ii,], edge.arrow.size = 2, edge.width = 4, edge.color = "black")
  dev.off()
  
  co_com <- layout_with_fr(n_com,
                       minx = rep(-Inf, vcount(n_com)),
                       miny = rep(-Inf, vcount(n_com)),
                       maxx = rep(Inf, vcount(n_com)),
                       maxy = rep(Inf, vcount(n_com)))
  filename <- paste0("com_topo_", sc,"_network.png")
  png(filename, width = 1000, height = 1000)
  plot(n_com, layout = co_com,
       edge.arrow.size = 2,
       edge.width = 4,
       edge.color = "black",
       vertex.label = NA)
  dev.off()
  
  
  Nodestats <- bind_rows(Nodestats,
                         tibble(node = V(n_com)$name,
                                degree_in =   degree(n_com, mode = "in"),
                                degree_out = degree(n_com, mode = "out"),
                                betweenness = betweenness(n_com),
                                closeness = closeness(n_com),
                                eigen_centrality = eigen_centrality(n_com)$vector,
                                page_rank = page_rank(n_com)$vector,
                                Community = sc))
  
  Graphstats <- bind_rows(Graphstats,
                          tibble(eigen_centrality = eigen_centrality(n_com)$value,
                                 page_rank = page_rank(n_com)$value,
                                 edge_density = edge_density(n_com),
                                 diameter = diameter(n_com),
                                 Community = sc))
}

#' Degree
Nodestats %>%
  mutate(Group = "Temperature sensitive") %>%
  mutate(Group = replace(Group, Community %in% paste0("R",7:12), "Membership sensitive")) %>%
  group_by(Community) %>%
  summarise(mean_degree_in = mean(degree_in),
            Group = unique(Group),
            .groups = "drop") %>%
  ggplot(aes(x = Group, y = mean_degree_in, col = Group)) + 
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  scale_color_manual(values = c("#d73027","#4575b4")) +
  
  theme_classic() +
  theme(axis.text.x = element_blank())


#' Betweenness
Nodestats %>%
  mutate(Group = "Temperature sensitive") %>%
  mutate(Group = replace(Group,
                         Community %in% paste0("R",7:12), "Membership sensitive")) %>%
  # group_by(Community) %>%
  # summarise(mean_degree_in = mean(degree_in),
  #           Group = unique(Group),
  #           .groups = "drop") %>%
  ggplot(aes(x = Group, y = betweenness, col = Group)) + 
  geom_point(position = position_jitterdodge()) +
  scale_color_manual(values = c("#d73027","#4575b4")) +
  
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_blank())


#' Closeness 
p1 <- Nodestats %>%
  mutate(Group = "Temperature sensitive") %>%
  mutate(Group = replace(Group,
                         Community %in% paste0("R",7:12), "Membership sensitive")) %>%
  # group_by(Community) %>%
  # summarise(mean_degree_in = mean(degree_in),
  #           Group = unique(Group),
  #           .groups = "drop") %>%
  ggplot(aes(x = Group, y = closeness, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), size = 3) +
  scale_color_manual(values = c("#d73027","#4575b4")) +
  guides(color = guide_legend(title = "Community type")) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
p1
ggsave("closeness.png", width = 6, height = 4)
ggsave("closeness.svg", width = 6, height = 4)

 #' Eigen
Nodestats %>%
  mutate(Group = "Temperature sensitive") %>%
  mutate(Group = replace(Group,
                         Community %in% paste0("R",7:12), "Membership sensitive")) %>%
  # group_by(Community) %>%
  # summarise(mean_degree_in = mean(degree_in),
  #           Group = unique(Group),
  #           .groups = "drop") %>%
  ggplot(aes(x = Group, y = eigen_centrality, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  scale_color_manual(values = c("#d73027","#4575b4")) +
  
  theme_classic() +
  theme(axis.text.x = element_blank())

#' Page rank
Nodestats %>%
  mutate(Group = "Temperature sensitive") %>%
  mutate(Group = replace(Group,
                         Community %in% paste0("R",7:12), "Membership sensitive")) %>%
  # group_by(Community) %>%
  # summarise(mean_degree_in = mean(degree_in),
  #           Group = unique(Group),
  #           .groups = "drop") %>%
  ggplot(aes(x = Group, y = page_rank, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  scale_color_manual(values = c("#d73027","#4575b4")) +
  
  theme_classic() +
  theme(axis.text.x = element_blank())






#' Eigen graph
Graphstats %>%
  mutate(Group = "Temperature sensitive") %>%
  mutate(Group = replace(Group,
                         Community %in% paste0("R",7:12), "Membership sensitive")) %>%
  ggplot(aes(x = Group, y = eigen_centrality, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  scale_color_manual(values = c("#d73027","#4575b4")) +
  
  theme_classic() +
  theme(axis.text.x = element_blank())

#' edge density
Graphstats %>%
  mutate(Group = "Temperature sensitive") %>%
  mutate(Group = replace(Group,
                         Community %in% paste0("R",7:12), "Membership sensitive")) %>%
  ggplot(aes(x = Group, y = edge_density, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  scale_color_manual(values = c("#d73027","#4575b4")) +
  
  theme_classic() +
  theme(axis.text.x = element_blank())

#' diameter
p1 <- Graphstats %>%
  mutate(Group = "Temperature sensitive") %>%
  mutate(Group = replace(Group,
                         Community %in% paste0("R",7:12), "Membership sensitive")) %>%
  ggplot(aes(x = Group, y = diameter, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), size = 3) +
  scale_color_manual(values = c("#d73027","#4575b4")) +
  guides(color = guide_legend(title = "Community type")) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
p1
ggsave("diameter.png", width = 6, height = 4)
ggsave("diameter.svg", width = 6, height = 4)
