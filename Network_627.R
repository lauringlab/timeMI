library(tidyverse)
library(associationsubgraphs)


mi_network_info <- readRDS("data/mi_network_info.rds")

pol_meaned_mi <- readRDS("data/MI_results_091122/h3n2_meaned_mi_noHA.rds") %>%
  left_join(mi_network_info, by = c("V1" = "index")) %>%
  left_join(mi_network_info, by = c("V2" = "index")) %>%
  select(a = "id.x", b = "id.y", strength = mip)
## Protein length cut-off
co <- 2232/2

pol_meaned_mi_co <- pol_meaned_mi %>%
  slice_max(strength, n = co)

net627_degree1 <- pol_meaned_mi_co %>%
  filter(a == "PB2_627" |
           b == "PB2_627")

deg1_ids <- unique(c(net627_degree1$a, net627_degree1$b))

net627_degree2 <- pol_meaned_mi_co %>%
  filter(a %in% deg1_ids |
           b %in% deg1_ids)

visualize_association_network(net627_degree1,
                             node_info = mi_network_info)

### With HA

pol_meaned_mi_HA <- readRDS("data/MI_results_091122/h3n2_meaned_mi.rds") %>%
  left_join(mi_network_info, by = c("V1" = "index")) %>%
  left_join(mi_network_info, by = c("V2" = "index")) %>%
  select(a = "id.x", b = "id.y", strength = mip)
## Protein length cut-off
coHA <- 2787/2

pol_meaned_mi_HA_co <- pol_meaned_mi_HA %>%
  slice_max(strength, n = coHA)

net627_degree1_HA <- pol_meaned_mi_HA_co %>%
  filter(a == "PB2_627" |
           b == "PB2_627")

deg1_ids_HA <- c(net627_degree1_HA$a, net627_degree1_HA$b)

visualize_association_network(net627_degree1,
                              node_info = mi_network_info)


# Z score 
mi_mean <- mean(pol_meaned_mi$strength)
mi_sd <- sd(pol_meaned_mi$strength)
pol_z <- pol_meaned_mi %>%
  mutate(strength_z = (strength - mi_mean)/mi_sd,
         id = paste(a,b, sep = "_"))

pol_z_top <- pol_z %>%
  filter(strength_z > 4)

net627_degree1z <- pol_z_top %>%
  filter(a == "PB2_627" |
           b == "PB2_627")

deg1_idsz <- unique(c(net627_degree1z$a, net627_degree1z$b))

