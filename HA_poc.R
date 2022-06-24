#This will explore HA to show proof of concept
library(tidyverse)
library(bio3d)


pdb_ha <- read.pdb("4o5n")
ha_ca.inds <- atom.select(pdb_ha, "calpha")
ha_dm <- dm(pdb_ha, ha_ca.inds)



pb2 <- readRDS("../H3N2_Cooccurrence/data/pb2_length.rds")
pb1 <- readRDS("../H3N2_Cooccurrence/data/pb1_length.rds")
pa <- readRDS("../H3N2_Cooccurrence/data/pa_length.rds")
ha <- readRDS("../H3N2_Cooccurrence/data/ha_length.rds")

apc <- readRDS("data/apc_hageog_imputed.rds") %>%
  filter(V1 > (pb2+pb1+pa) & V2 > (pb2+pb1+pa)) %>%
  mutate(V1 = V1-(pb2+pb1+pa),
         V2 = V2-(pb2+pb1+pa))
meanMI <- mean(apc$mi)

apc_unweighted <- readRDS("data/apc_ha_imputed.rds") %>%
  filter(V1 > (pb2+pb1+pa) & V2 > (pb2+pb1+pa)) %>%
  mutate(V1 = V1-(pb2+pb1+pa),
         V2 = V2-(pb2+pb1+pa))
meanMI_unweighted <- mean(apc_unweighted$mi)

MIp <- apc %>%
  mutate(apc = (v1_meanMIgeog*v2_meanMIgeog)/meanMI,
         mip = mi - apc) %>%
  mutate(Group = paste("[", V1, ";", V2, "]", sep = "")) %>%
  select(Group,mi, mip, apc, V1, V2) %>%
  pivot_longer(cols = c(V1, V2), names_to = "x_type", values_to = "x_val") %>%
  mutate(sector = "HA")

MIp_unweighted <- apc_unweighted %>%
  mutate(apc = (v1_meanMI*v2_meanMI)/meanMI_unweighted,
         mip = mi - apc) %>%
  mutate(Group = paste("[", V1, ";", V2, "]", sep = "")) %>%
  select(Group,mi, mip, apc, V1, V2) %>%
  pivot_longer(cols = c(V1, V2), names_to = "x_type", values_to = "x_val") %>%
  mutate(sector = "HA")


get_res_indices <- function(atom){
  atom %>%
    filter(elety == "CA") %>%
    select(resid, chain, resno) %>%
    mutate(mat_index = 1:length(.$resid)) %>%
    mutate(resno = case_when(
      chain == "B" ~ as.numeric(resno + 343),
      TRUE ~ as.numeric(resno)
    ))
}

ha_res_indices <- get_res_indices(pdb_ha$atom)
get_distances <- function(mips, res_indices, pdb_dist) {
  dict <- mips %>%
    full_join(res_indices, by = c("x_val" = "resno")) %>%
    filter(!is.na(mip))
  
  dist <- dict %>%
    select(Group, chain, mat_index, x_type) %>%
    filter(!is.na(mat_index)) %>%
    pivot_wider(names_from = x_type, values_from = mat_index) %>%
    filter(!is.na(V1) & !is.na(V2)) %>%
    group_by(Group) %>%
    mutate(minx = min(V1, V2),
           maxx = max(V1, V2)) %>%
    mutate(distance = pdb_dist[minx, maxx])
  
  distances <- dict %>%
    full_join(dist, by = "Group") %>%
    ungroup() %>%
    select(Group, mi, mip, apc, distance) %>%
    filter(!is.na(Group)) %>%
    unique() %>%
    filter(!is.na(distance))
  
  return(distances)
}

ha_distances <- get_distances(MIp, ha_res_indices, ha_dm)
ha_distances_unweighted <- get_distances(MIp_unweighted, ha_res_indices, ha_dm)


nums5 <- seq(5, 300, 5)

get_mean_dists <- function(nums, distances, order){
  order = enquo(order)
  mean_dists <- c()
  for (i in nums){
    mean_dists = c(mean_dists,
                   distances %>%
                     slice_max(n = i, order_by = !!order) %>%
                     pull(distance) %>%
                     mean())
  }
  return(mean_dists)
}

ha_mean_dist_mi <- get_mean_dists(nums5, ha_distances, order = mi)
ha_mean_dist_mip <- get_mean_dists(nums5, ha_distances, order = mip)

ha_mean_dist_mi_unweighted <- get_mean_dists(nums5, ha_distances_unweighted, order = mi)
ha_mean_dist_mip_unweighted <- get_mean_dists(nums5, ha_distances_unweighted, order = mip)

dist_compare_ha <- tibble(
  num = nums5,
  unweighted_mi = ha_mean_dist_mi_unweighted,
  unweighted_mip = ha_mean_dist_mip_unweighted,
  weighted_mi = ha_mean_dist_mi,
  weighted_mip = ha_mean_dist_mip
) %>%
  pivot_longer(cols = !num, names_to = "analysis", values_to = "distance")

ggplot(dist_compare_ha, aes(x = num, y = distance, group = analysis, color = analysis)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5)

##################################Above 20######################################
get_distances_above20 <- function(mips, res_indices, pdb_dist) {
  dict <- mips %>%
    full_join(res_indices, by = c("x_val" = "resno")) %>%
    filter(!is.na(mip))
  
  dist <- dict %>%
    select(Group, chain, mat_index, x_type) %>%
    filter(!is.na(mat_index)) %>%
    pivot_wider(names_from = x_type, values_from = mat_index) %>%
    filter(!is.na(V1) & !is.na(V2)) %>%
    filter(abs(V2-V1) > 20) %>%
    group_by(Group) %>%
    mutate(minx = min(V1, V2),
           maxx = max(V1, V2)) %>%
    mutate(distance = pdb_dist[minx, maxx])
  
  distances <- dict %>%
    full_join(dist, by = "Group") %>%
    select(Group, mi, mip, apc, distance) %>%
    unique() %>%
    filter(!is.na(distance))
  
  return(distances)
}

ha_distances_above20 <- get_distances_above20(MIp, ha_res_indices, ha_dm)
ha_distances_unweighted_above20 <- get_distances_above20(MIp_unweighted, ha_res_indices, ha_dm)


ha_mean_dist_mi_above20 <- get_mean_dists(nums5, ha_distances_above20, order = mi)
ha_mean_dist_mip_above20 <- get_mean_dists(nums5, ha_distances_above20, order = mip)

ha_mean_dist_mi_unweighted_above20 <- get_mean_dists(nums5, ha_distances_unweighted_above20, order = mi)
ha_mean_dist_mip_unweighted_above20 <- get_mean_dists(nums5, ha_distances_unweighted_above20, order = mip)

dist_compare_ha_above20 <- tibble(
  num = nums5,
  unweighted_mi = ha_mean_dist_mi_unweighted_above20,
  unweighted_mip = ha_mean_dist_mip_unweighted_above20,
  weighted_mi = ha_mean_dist_mi_above20,
  weighted_mip = ha_mean_dist_mip_above20
) %>%
  pivot_longer(cols = !num, names_to = "analysis", values_to = "distance")

ggplot(dist_compare_ha_above20, aes(x = num, y = distance, group = analysis, color = analysis)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5)




