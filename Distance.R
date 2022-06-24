library(bio3d)
library(tidyverse)

pb2 <- readRDS("../H3N2_Cooccurrence/data/pb2_length.rds")
pb1 <- readRDS("../H3N2_Cooccurrence/data/pb1_length.rds")
pa <- readRDS("../H3N2_Cooccurrence/data/pa_length.rds")

pdb_prei <- read.pdb("4wsb")
pdb_pcs <- read.pdb("6rr7")
pdb_elo <- read.pdb("6t0v")
pdb_ter <- read.pdb("6szu")

# prei_ca.inds <- atom.select(pdb_prei, "calpha")
# prei_dm <- dm(pdb_prei, prei_ca.inds)  
# saveRDS(prei_dm, "4wsb_dm.rda")
# 
# elo_ca.inds <- atom.select(pdb_elo, "calpha")
# elo_dm <- dm(pdb_elo, elo_ca.inds)  
# saveRDS(elo_dm, "6t0v_dm.rda")
#
# ter_ca.inds <- atom.select(pdb_ter, "calpha")
# ter_dm <- dm(pdb_ter, ter_ca.inds)
# saveRDS(ter_dm, "6szu_dm.rda")

prei_dist <- readRDS("4wsb_dm.rda")
pcs_dist <- readRDS("../H3N2_Cooccurrence/data/6rr7_dm.rda")
elo_dist <- readRDS("6t0v_dm.rda")
ter_dist <- readRDS("6szu_dm.rda")
#Large C-terminal chunk of PB2 is unmodeled in the termination structure

apc <- readRDS("data/apc_hageog_imputed.rds") %>%
  filter(V1 < (pb2+pb1+pa) & V2 < (pb2+pb1+pa))
meanMI <- mean(apc$mi)

apc_unweighted <- readRDS("data/apc_ha_imputed.rds") %>%
  filter(V1 < (pb2+pb1+pa) & V2 < (pb2+pb1+pa))
meanMI_unweighted <- mean(apc_unweighted$mi)

MIp <- apc %>%
  mutate(apc = (v1_meanMIgeog*v2_meanMIgeog)/meanMI,
         mip = mi - apc) %>%
  mutate(Group = paste("[", V1, ";", V2, "]", sep = "")) %>%
  select(Group,mi, mip, apc, V1, V2) %>%
  pivot_longer(cols = c(V1, V2), names_to = "x_type", values_to = "x_val") %>%
  mutate(sector = case_when(
    x_val <= pb2 ~ "PB2",
    x_val <= pb2+pb1 ~ "PB1",
    TRUE ~ "PA"
  )) %>%
  mutate(x_val = case_when(
    sector == "PB1" ~ as.numeric(x_val - pb2),
    sector == "PA" ~ as.numeric(x_val - (pb2+pb1)),
    TRUE ~ as.numeric(x_val)
  ))

MIp_unweighted <- apc_unweighted %>%
  mutate(apc = (v1_meanMI*v2_meanMI)/meanMI_unweighted,
         mip = mi - apc) %>%
  mutate(Group = paste("[", V1, ";", V2, "]", sep = "")) %>%
  select(Group,mi, mip, apc, V1, V2) %>%
  pivot_longer(cols = c(V1, V2), names_to = "x_type", values_to = "x_val") %>%
  mutate(sector = case_when(
    x_val <= pb2 ~ "PB2",
    x_val <= pb2+pb1 ~ "PB1",
    TRUE ~ "PA"
  )) %>%
  mutate(x_val = case_when(
    sector == "PB1" ~ as.numeric(x_val - pb2),
    sector == "PA" ~ as.numeric(x_val - (pb2+pb1)),
    TRUE ~ as.numeric(x_val)
  ))

get_res_indices <- function(atom){
  atom %>%
    filter(elety == "CA") %>%
    select(resid, chain, resno) %>%
    mutate(mat_index = 1:length(.$resid)) %>%
    mutate(sector = case_when(
      chain == "A" ~ "PA",
      chain == "B" ~ "PB1",
      TRUE ~ "PB2"
    ))
}

prei_res_indices <- get_res_indices(pdb_prei$atom)
pcs_res_indices <- get_res_indices(pdb_pcs$atom)
elo_res_indices <- get_res_indices(pdb_elo$atom)
ter_res_indices <- get_res_indices(pdb_ter$atom)

get_distances <- function(mips, res_indices, pdb_dist) {
  dict <- mips %>%
    full_join(res_indices, by = c("sector", "x_val" = "resno")) %>%
    filter(!is.na(mip))
  
  dist <- dict %>%
    select(Group, mat_index, x_type) %>%
    pivot_wider(names_from = x_type, values_from = mat_index) %>%
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

prei_distances <- get_distances(MIp, prei_res_indices, prei_dist)
pcs_distances <- get_distances(MIp, pcs_res_indices, pcs_dist)
elo_distances <- get_distances(MIp, elo_res_indices, elo_dist)
ter_distances <- get_distances(MIp, ter_res_indices, ter_dist)

prei_distances_unweighted <- get_distances(MIp_unweighted, prei_res_indices, prei_dist)
pcs_distances_unweighted <- get_distances(MIp_unweighted, pcs_res_indices, pcs_dist)
elo_distances_unweighted <- get_distances(MIp_unweighted, elo_res_indices, elo_dist)
ter_distances_unweighted <- get_distances(MIp_unweighted, ter_res_indices, ter_dist)

# ggplot(ter_distances, aes(x = distance, y = mip)) +
#   geom_point(alpha = 0.3)

#####################Plot of mean distance by top residues######################

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

prei_mean_dist_mi <- get_mean_dists(nums5, prei_distances, order = mi)
prei_mean_dist_mip <- get_mean_dists(nums5, prei_distances, order = mip)
pcs_mean_dist_mi <- get_mean_dists(nums5, pcs_distances, order = mi)
pcs_mean_dist_mip <- get_mean_dists(nums5, pcs_distances, order = mip)
elo_mean_dist_mi <- get_mean_dists(nums5, elo_distances, order = mi)
elo_mean_dist_mip <- get_mean_dists(nums5, elo_distances, order = mip)
ter_mean_dist_mi <- get_mean_dists(nums5, ter_distances, order = mi)
ter_mean_dist_mip <- get_mean_dists(nums5, ter_distances, order = mip)

prei_mean_dist_mi_unweighted <- get_mean_dists(nums5, prei_distances_unweighted, order = mi)
prei_mean_dist_mip_unweighted <- get_mean_dists(nums5, prei_distances_unweighted, order = mip)
pcs_mean_dist_mi_unweighted <- get_mean_dists(nums5, pcs_distances_unweighted, order = mi)
pcs_mean_dist_mip_unweighted <- get_mean_dists(nums5, pcs_distances_unweighted, order = mip)
elo_mean_dist_mi_unweighted <- get_mean_dists(nums5, elo_distances_unweighted, order = mi)
elo_mean_dist_mip_unweighted <- get_mean_dists(nums5, elo_distances_unweighted, order = mip)
ter_mean_dist_mi_unweighted <- get_mean_dists(nums5, ter_distances_unweighted, order = mi)
ter_mean_dist_mip_unweighted <- get_mean_dists(nums5, ter_distances_unweighted, order = mip)



####################Now plot against non-weighted MIp and MI####################
dist_compare_prei <- tibble(
  num = nums5,
  unweighted_mi = prei_mean_dist_mi_unweighted,
  unweighted_mip = prei_mean_dist_mip_unweighted,
  weighted_mi = prei_mean_dist_mi,
  weighted_mip = prei_mean_dist_mip
) %>%
  pivot_longer(cols = !num, names_to = "analysis", values_to = "distance")

dist_compare_pcs <- tibble(
  num = nums5,
  unweighted_mi = pcs_mean_dist_mi_unweighted,
  unweighted_mip = pcs_mean_dist_mip_unweighted,
  weighted_mi = pcs_mean_dist_mi,
  weighted_mip = pcs_mean_dist_mip
) %>%
  pivot_longer(cols = !num, names_to = "analysis", values_to = "distance")

dist_compare_elo <- tibble(
  num = nums5,
  unweighted_mi = elo_mean_dist_mi_unweighted,
  unweighted_mip = elo_mean_dist_mip_unweighted,
  weighted_mi = elo_mean_dist_mi,
  weighted_mip = elo_mean_dist_mip
) %>%
  pivot_longer(cols = !num, names_to = "analysis", values_to = "distance")

dist_compare_ter <- tibble(
  num = nums5,
  unweighted_mi = ter_mean_dist_mi_unweighted,
  unweighted_mip = ter_mean_dist_mip_unweighted,
  weighted_mi = ter_mean_dist_mi,
  weighted_mip = ter_mean_dist_mip
) %>%
  pivot_longer(cols = !num, names_to = "analysis", values_to = "distance")

dist_compare_all <- bind_rows(list(
  "Pre-Initiaion" = dist_compare_prei,
  "Post Cap-snatching" = dist_compare_pcs,
  "Elongation" = dist_compare_elo,
  "Termination" = dist_compare_ter
), .id = "Structure")

saveRDS(dist_compare_all, "dist_compare_all.rds")

ggplot(dist_compare_all, aes(x = num, y = distance, group = analysis, color = analysis)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  facet_wrap(~Structure)


##########Do analysis with only residues further than 20 residues apart#########

get_distances_above20 <- function(mips, res_indices, pdb_dist) {
  dict <- mips %>%
    full_join(res_indices, by = c("sector", "x_val" = "resno")) %>%
    filter(!is.na(mip))
  
  dist <- dict %>%
    select(Group, mat_index, x_type) %>%
    pivot_wider(names_from = x_type, values_from = mat_index) %>%
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

prei_distances_above20 <- get_distances_above20(MIp, prei_res_indices, prei_dist)
pcs_distances_above20 <- get_distances_above20(MIp, pcs_res_indices, pcs_dist)
elo_distances_above20 <- get_distances_above20(MIp, elo_res_indices, elo_dist)
ter_distances_above20 <- get_distances_above20(MIp, ter_res_indices, ter_dist)

prei_distances_above20_unweighted <- get_distances_above20(MIp_unweighted, prei_res_indices, prei_dist)
pcs_distances_above20_unweighted <- get_distances_above20(MIp_unweighted, pcs_res_indices, pcs_dist)
elo_distances_above20_unweighted <- get_distances_above20(MIp_unweighted, elo_res_indices, elo_dist)
ter_distances_above20_unweighted <- get_distances_above20(MIp_unweighted, ter_res_indices, ter_dist)

prei_above20_mean_dist_mi <- get_mean_dists(nums5, prei_distances_above20, order = mi)
prei_above20_mean_dist_mip <- get_mean_dists(nums5, prei_distances_above20, order = mip)
pcs_above20_mean_dist_mi <- get_mean_dists(nums5, pcs_distances_above20, order = mi)
pcs_above20_mean_dist_mip <- get_mean_dists(nums5, pcs_distances_above20, order = mip)
elo_above20_mean_dist_mi <- get_mean_dists(nums5, elo_distances_above20, order = mi)
elo_above20_mean_dist_mip <- get_mean_dists(nums5, elo_distances_above20, order = mip)
ter_above20_mean_dist_mi <- get_mean_dists(nums5, ter_distances_above20, order = mi)
ter_above20_mean_dist_mip <- get_mean_dists(nums5, ter_distances_above20, order = mip)

prei_above20_mean_dist_mi_unweighted <- get_mean_dists(nums5, prei_distances_above20_unweighted, order = mi)
prei_above20_mean_dist_mip_unweighted <- get_mean_dists(nums5, prei_distances_above20_unweighted, order = mip)
pcs_above20_mean_dist_mi_unweighted <- get_mean_dists(nums5, pcs_distances_above20_unweighted, order = mi)
pcs_above20_mean_dist_mip_unweighted <- get_mean_dists(nums5, pcs_distances_above20_unweighted, order = mip)
elo_above20_mean_dist_mi_unweighted <- get_mean_dists(nums5, elo_distances_above20_unweighted, order = mi)
elo_above20_mean_dist_mip_unweighted <- get_mean_dists(nums5, elo_distances_above20_unweighted, order = mip)
ter_above20_mean_dist_mi_unweighted <- get_mean_dists(nums5, ter_distances_above20_unweighted, order = mi)
ter_above20_mean_dist_mip_unweighted <- get_mean_dists(nums5, ter_distances_above20_unweighted, order = mip)

dist_compare_prei_above20 <- tibble(
  num = nums5,
  unweighted_mi = prei_above20_mean_dist_mi_unweighted,
  unweighted_mip = prei_above20_mean_dist_mip_unweighted,
  weighted_mi = prei_above20_mean_dist_mi,
  weighted_mip = prei_above20_mean_dist_mip
) %>%
  pivot_longer(cols = !num, names_to = "analysis", values_to = "distance")

dist_compare_pcs_above20 <- tibble(
  num = nums5,
  unweighted_mi = pcs_above20_mean_dist_mi_unweighted,
  unweighted_mip = pcs_above20_mean_dist_mip_unweighted,
  weighted_mi = pcs_above20_mean_dist_mi,
  weighted_mip = pcs_above20_mean_dist_mip
) %>%
  pivot_longer(cols = !num, names_to = "analysis", values_to = "distance")

dist_compare_elo_above20 <- tibble(
  num = nums5,
  unweighted_mi = elo_above20_mean_dist_mi_unweighted,
  unweighted_mip = elo_above20_mean_dist_mip_unweighted,
  weighted_mi = elo_above20_mean_dist_mi,
  weighted_mip = elo_above20_mean_dist_mip
) %>%
  pivot_longer(cols = !num, names_to = "analysis", values_to = "distance")

dist_compare_ter_above20 <- tibble(
  num = nums5,
  unweighted_mi = ter_above20_mean_dist_mi_unweighted,
  unweighted_mip = ter_above20_mean_dist_mip_unweighted,
  weighted_mi = ter_above20_mean_dist_mi,
  weighted_mip = ter_above20_mean_dist_mip
) %>%
  pivot_longer(cols = !num, names_to = "analysis", values_to = "distance")

dist_compare_all_above20 <- bind_rows(list(
  "Pre-Initiaion" = dist_compare_prei_above20,
  "Post Cap-snatching" = dist_compare_pcs_above20,
  "Elongation" = dist_compare_elo_above20,
  "Termination" = dist_compare_ter_above20
), .id = "Structure")

saveRDS(dist_compare_all_above20, "dist_compare_all_above20.rds")


ggplot(dist_compare_all_above20, aes(x = num, y = distance, group = analysis, color = analysis)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  facet_wrap(~Structure)





##########################De-bugging residues indices###########################

seq_times <- readRDS("/Users/saraharcos/Desktop/Lauring Lab/MI_Networks_App/data/seq_times.rda")

seq1 <- seq_times[1,]$x %>%
  str_split(pattern = "") %>%
  unlist()

codes <- tibble(
  res = c("A", "R", "N", "D", "C",
          "E", "Q", "G", "H", "I",
          "L", "K", "M", "F", "P",
          "S", "T", "W", "Y", "V",
          "-", "X"),
  res_full = c("ALA", "ARG", "ASN", "ASP", "CYS",
               "GLU", "GLN", "GLY", "HIS", "ILE",
               "LEU", "LYS", "MET", "PHE", "PRO",
               "SER", "THR", "TRP", "TYR", "VAL",
               "GAP", "UNK")
)

debug <- tibble(
  sector = c(rep("PB2", pb2), rep("PB1", pb1), rep("PA", pa)),
  resno = c(1:pb2, 1:pb1, 1:pa),
  res = seq1[1:(pb2+pb1+pa)]
) %>%
  left_join(codes, by = "res")

test <- res_indices %>%
  full_join(debug, by = c("resno", "sector"))

problems <- test %>%
  filter(resid != res_full)

p <- seq_times %>%
  mutate(x2 = str_sub(x, start = 692, end = 692))












