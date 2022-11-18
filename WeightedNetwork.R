library(tidyverse)
library(associationsubgraphs)


mi_network_info <- tibble(
  index = 1:2798,
  gene = c(rep("PB2", 759), rep("PB1", 757), rep("PA", 716), rep("HA", 566)),
  id = c(1:759, 1:757, 1:716, 1:566)) %>%
  mutate(id = paste(gene, id, sep = "_"),
         color = case_when(
    str_detect(id, "PA") ~ "#00A14B",
    str_detect(id, "PB1") ~ "#08519c",
    str_detect(id, "PB2") ~ "#e41a1c",
    TRUE ~ "#fd8d3c"
  )) %>%
  select(-gene)

saveRDS(mi_network_info, "data/mi_network_info.rds")

mi <- readRDS("data/MI_results_091122/h3n2_meaned_mi.rds") %>%
  select(V1, V2, mip) %>%
  left_join(mi_network_info, by = c("V1" = "index")) %>%
  left_join(mi_network_info, by = c("V2" = "index")) %>%
  select(a = "id.x", b = "id.y", strength = mip)

mi_mean <- mean(mi$strength)
mi_sd <- sd(mi$strength)

mi_z <- mi %>%
  mutate(strength = (strength - mi_mean)/mi_sd) %>%
  filter(strength > 4)

mi_noHA<- readRDS("data/MI_results_091122/h3n2_meaned_mi_noHA.rds") %>%
  left_join(mi_network_info, by = c("V1" = "index")) %>%
  left_join(mi_network_info, by = c("V2" = "index")) %>%
  select(a = "id.x", b = "id.y", strength = mip)

mi_noHA_mean <- mean(mi_noHA$strength)
mi_noHA_sd <- sd(mi_noHA$strength)

mi_noHA_z <- mi_noHA %>%
  mutate(strength = (strength - mi_noHA_mean)/mi_noHA_sd) %>%
  filter(strength > 4)





### Within between pair tallys ######

mi_noHA_tallys <- mi_noHA_z %>%
  mutate(grouping = paste(a,b, sep = ",")) %>%
  mutate(grouping = str_remove_all(grouping, pattern = "_[0-9]*")) %>%
  count(grouping) %>%
  arrange(n) %>% 
  mutate(grouping = factor(grouping, levels = unique(grouping)))

ggplot(mi_noHA_tallys, aes(x = grouping, y = n)) +
  geom_col(alpha = 0.7, fill = "grey") +
  theme_classic() +
  labs(title = "wMI pairs within and between polymerase proteins",
       x = "Group", y = "Count") +
  theme(text = element_text(size = 12),
        legend.position = "none")

ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/mip_pairs.pdf",
       width = 5, height = 3)

mi_tallys <- mi_z %>%
  mutate(group1 = paste(a,b, sep = ",")) %>%
  mutate(group1 = str_remove_all(group1, pattern = "_[0-9]*")) %>%
  mutate(grouping = case_when(
    group1 %in% c("PB1,HA", "PB2,HA", "PA,HA") ~ "Pol-HA",
    group1 == "HA,HA" ~ "HA-HA",
    TRUE ~ "Pol-Pol"
  )) %>%
  count(grouping) %>%
  arrange(n) %>% 
  mutate(grouping = factor(grouping, levels = unique(grouping)))

ggplot(mi_tallys, aes(x = grouping, y = n)) +
  geom_col(alpha = 0.7, fill = "grey") +
  theme_classic() +
  labs(title = "wMI pairs within and between the polymerase and HA",
       x = "Group", y = "Count") +
  theme(text = element_text(size = 12),
        legend.position = "none")


ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/pol_HA_pairs.pdf",
       width = 3, height = 3)

#Calculate min_max_rule with full data first, then set default_step but only
#Use top pairs

visualize_subgraph_structure(
  association_pairs = mi_noHA,
  node_info = mi_network_info %>% select(-index),
  trim_subgraph_results = FALSE,
  default_step = "min_max_rule"
)

visualize_subgraph_structure(
  association_pairs = mi,
  node_info = mi_network_info %>% select(-index),
  trim_subgraph_results = FALSE,
  default_step = "min_max_rule"
)


##Large dataset for presentation with hairball


mi_noHA_big<- readRDS("data/MI_results_091122/h3n2_meaned_mi_noHA.rds") %>%
  left_join(mi_network_info, by = c("V1" = "index")) %>%
  left_join(mi_network_info, by = c("V2" = "index")) %>%
  select(a = "id.x", b = "id.y", strength = mip) %>%
  slice_max(strength, n = 5000)

visualize_subgraph_structure(
  association_pairs = mi_noHA_big,
  node_info = mi_network_info %>% select(-index),
  trim_subgraph_results = FALSE,
  default_step = "min_max_rule"
)

visualize_association_network(association_pairs = mi_noHA_big,
                              node_info = mi_network_info %>% select(-index))

#Recalculate MIp without HA

pb2 <- readRDS("/Users/saraharcos/Desktop/Lauring Lab/H3N2_Cooccurrence/data/pb2_length.rds")
pb1 <- readRDS("/Users/saraharcos/Desktop/Lauring Lab/H3N2_Cooccurrence/data/pb1_length.rds")
pa <- readRDS("/Users/saraharcos/Desktop/Lauring Lab/H3N2_Cooccurrence/data/pa_length.rds")

mi_noHA <- readRDS("data/h3n2_meaned_mi.rds") %>%
  filter(V1 < (pb2+pb1+pa)) %>%
  filter(V2 < (pb2+pb1+pa))



meanMI <- mean(mi_noHA$mi)

indices <- unique(c(mi_noHA$V1, mi_noHA$V2))

mean_MIs <- numeric(length(indices))
for (i in indices){
  mean_MIs[i] <- mi_noHA[mi_noHA$V1 == i | mi_noHA$V2 == i,] %>% pull(mi) %>% mean()
}

v1_meanMI <- numeric(nrow(mi_noHA))
v2_meanMI <- numeric(nrow(mi_noHA))
for (i in 1:nrow(mi_noHA)) {
  v1_meanMI[i] <- mean_MIs[mi_noHA$V1[i]]
  v2_meanMI[i] <- mean_MIs[mi_noHA$V2[i]]
}

mi_noHA$v1_meanMI <- v1_meanMI
mi_noHA$v2_meanMI <- v2_meanMI

mip_noHA <- mi_noHA %>% mutate(apc = (v1_meanMI*v2_meanMI)/meanMI,
                             mip = mi - apc) %>%
         mutate(Group = paste("[", V1, ";", V2, "]", sep = ""))

saveRDS(mip_noHA, "data/h3n2_meaned_mip_noHA.rds")




