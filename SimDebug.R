#Testing code with simulated data
#Simulated data generated using iqtree2 --alisim alignment_ph -t RANDOM{yh/7000} -m LG --length 2800

library(tidyverse)
library(purrr)
library(data.table)
library(foreach)
library(doSNOW)
library(readxl)
library(lubridate)

spike_times <- readRDS("data/spike_times.rds") %>%
  arrange(date) %>%
  mutate(group = date)
seq_times <- readRDS("data/seq_times.rda")


spike_weights <- spike_times %>%
  select(date) %>%
  filter(!is.na(date)) %>%
  count(date) %>%
  mutate(freq = n/sum(n)) %>%
  mutate(groupnum = 1:12,
         group = date)

sampling_weights <- seq_times %>%
  select(date) %>%
  filter(!is.na(date)) %>%
  count(date) %>%
  mutate(sfreq = n/sum(n)) %>%
  top_n(n = 12, wt = date) %>%
  mutate(groupnum = 1:12) %>%
  select(groupnum, sfreq) %>%
  full_join(spike_weights, by = "groupnum")


run_sim <- function(phy){
  
  mat <- format_phy(phy)
  mi <- calculate_mi(mat)
  phy_s <- sample_phy(phy, sampling_weights)
  mat_s <- format_phy(phy_s)
  mi_s <- calculate_mi(mat_s)
  mi_sw <- calculate_weighted_mi(phy_s,mat_s, spike_weights)
  # saveRDS(list("mi" = mi, "mi_s" = mi_s, "mi_sw" = mi_sw),
  #         paste("data/sim_mi/", filename, ".rds", sep = ""))
  return(list("mi" = mi, "mi_s" = mi_s, "mi_sw" = mi_sw))
}

format_phy <- function(phy){
  sim_mat <- c(str_split(phy$x, "", simplify = TRUE))
  dim(sim_mat) <- c(1083,1274)
  rownames(sim_mat) = phy$rowname
  return(sim_mat)
}


calculate_mi <- function(mat){
  entropies <- vector(length = ncol(mat))
  for(i in 1:ncol(mat)){
    entropies[i] = get_ent_unweighted(i, mat)
  }
  
  non_zero_e <- which(entropies!=0)
  
  joint_cols <- t(combn(non_zero_e, 2)) %>%
    as.data.frame()
  
  joint_ents <- vector(length = nrow(joint_cols))
  for(i in 1:nrow(joint_cols)){
    joint_ents[i] = get_jointent_unweighted(joint_cols[i, "V1"], joint_cols[i, "V2"], mat)
  }
  
  return(format_mi(entropies, joint_ents, joint_cols))
}

format_mi <- function(entropies, joint_entropies, joint_cols){
  joint_cols$joint_entropy <- joint_entropies
  non_zero_e <- which(entropies!=0)
  
  v1_entropy <- numeric(nrow(joint_cols))
  v2_entropy <- numeric(nrow(joint_cols))
  for (i in 1:nrow(joint_cols)) {
    v1_entropy[i] <- entropies[joint_cols$V1[i]]
    v2_entropy[i] <- entropies[joint_cols$V2[i]]
  }
  
  joint_cols$v1_entropy <- v1_entropy
  joint_cols$v2_entropy <- v2_entropy
  
  joint_cols <- joint_cols %>%
    mutate(mi = v1_entropy + v2_entropy - joint_entropy)
  
  meanMI <- mean(joint_cols$mi)
  
  mean_MIs <- numeric(length(non_zero_e))
  for (i in non_zero_e){
    mean_MIs[i] <- joint_cols[joint_cols$V1 == i | joint_cols$V2 == i,] %>% pull(mi) %>% mean()
  }
  
  v1_meanMI <- numeric(nrow(joint_cols))
  v2_meanMI <- numeric(nrow(joint_cols))
  for (i in 1:nrow(joint_cols)) {
    v1_meanMI[i] <- mean_MIs[joint_cols$V1[i]]
    v2_meanMI[i] <- mean_MIs[joint_cols$V2[i]]
  }
  
  joint_cols$v1_meanMI <- v1_meanMI
  joint_cols$v2_meanMI <- v2_meanMI
  
  return(joint_cols %>% mutate(apc = (v1_meanMI*v2_meanMI)/meanMI,
                               mip = mi - apc) %>%
           mutate(Group = paste("[", V1, ";", V2, "]", sep = "")))
}

get_ent_unweighted <- function(residue, mat){
  t <- table(mat[,residue]) %>%
    prop.table()
  
  return(sum(-t*log(t)))
}

get_jointent_unweighted <- function(residue_a,residue_b, mat){
  jointcol <- as.matrix(paste(mat[,residue_a],mat[,residue_b],sep=""))
  
  t <- table(jointcol) %>%
    prop.table()
  
  return(sum(-t*log(t)))
}

sample_phy <- function(phy, weights){
  sampled_data <- tibble(
    x = phy$x,
    group = phy$group
  ) %>%
    left_join(weights, by = "group") %>%
    slice_sample(n = 1083, replace = TRUE, weight_by = sfreq) %>%
    mutate(rowname = c(1:1083))
}

calculate_weighted_mi <- function(df, mat, weights){
  inverse_weights <- weights %>%
    mutate(inverse_weight = 1/n) %>%
    mutate(inverse_weight = inverse_weight/sum(inverse_weight)) %>%
    select(-n)
  
  groups <- df %>%
    select(name, group) %>%
    nest(names = name) %>%
    full_join(inverse_weights, by = "group")
  
  entropies <- get_ent_weighted(mat, groups, inverse_weights)
  
  non_zero_e <- which(entropies!=0)
  
  joint_cols <- t(combn(non_zero_e, 2)) %>%
    as.data.frame()
  
  joint_ents <- vector(length = nrow(joint_cols))
  for(i in 1:nrow(joint_cols)) {
    joint_ents[i] <- get_weighted_group_jointent(joint_cols[i, "V1"], joint_cols[i, "V2"], mat, groups)
  }
  
  entropies <- unlist(entropies)
  return(format_mi(entropies, joint_ents, joint_cols))
}

get_weighted_group_freqs <- function(grp, mat, groups, inverse_weights){
  group_mat <- mat[unlist(groups[groups$group == grp,]$names),]
  
  group_weight <- inverse_weights %>%
    filter(group == grp) %>%
    pull(inverse_weight)
  
  if(is.null(dim(group_mat))){
    group_freq <- lapply(group_mat, table)
  }
  else {
    group_count <- apply(group_mat, 2, table)
    group_freq <- lapply(group_count, prop.table)
  }
  return(lapply(group_freq, function(x) as.list(x*group_weight)))
  
}

get_ent_weighted <- function(mat, groups, inverse_weights) {
  groups_vec <- unique(groups$group)
  
  weighted_group_freqs <- vector(mode = "list", length(groups_vec))
  weighted_group_freqs <- lapply(groups_vec, get_weighted_group_freqs, mat = mat,
                                 groups = groups,
                                 inverse_weights = inverse_weights)
  
  weighted_group_freqs_t <- vector(mode = "list", length(weighted_group_freqs[[2]]))
  weighted_group_freqs_t <- purrr::transpose(weighted_group_freqs)
  
  weighted_group_freqs_all <- vector(mode = "list", length(weighted_group_freqs_t))
  weighted_group_freqs_all<- lapply(weighted_group_freqs_t, rbindlist, fill = TRUE)
  
  weighted_group_ent <- vector(mode = "list", length(weighted_group_freqs_all))
  weighted_group_ent <- lapply( weighted_group_freqs_all, get_ent_helper)
  
  return(weighted_group_ent)
}

get_ent_helper <- function(pos){
  f <- colSums(pos, na.rm = TRUE)
  return(sum(-f*log(f)))
}

get_weighted_group_jointfreqs <- function(grp, mat, groups){
  group_mat <- mat[unlist(groups[groups$group == grp,]$names),]
  
  group_freq  <- group_mat %>% table() %>% prop.table()
  return(as.list(group_freq*groups[groups$group == grp,]$inverse_weight))
}

get_weighted_group_jointent <- function(colA, colB, mat, groups){
  joint_mat <- as.matrix(paste(mat[,colA],mat[,colB],sep=""))
  rownames(joint_mat) <- rownames(mat)
  
  groups_vec <- unique(groups$group)
  
  weighted_group_jointfreqs <- vector(mode = "list", length(groups_vec))
  for (i in 1:length(groups_vec)) {
    weighted_group_jointfreqs[[i]] = get_weighted_group_jointfreqs(groups_vec[i], joint_mat, groups)
  }
  
  return(rbindlist(weighted_group_jointfreqs,fill = TRUE) %>% get_ent_helper())
}


test <- run_sim(spike_times)

mip_compare <- test$mi %>%
  inner_join(test$mi_s, by = "Group") %>%
  inner_join(test$mi_sw, by = "Group")

cor(mip_compare$mip.x, mip_compare$mip.y)
cor(mip_compare$mip.x, mip_compare$mip)

ggplot(mip_compare, aes(x = mip.x, y = mip)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method='lm')



mat <- format_phy(spike_times)
mi <- calculate_mi(mat)
phy_s <- sample_phy(spike_times, sampling_weights)
mat_s <- format_phy(phy_s)
mi_s <- calculate_mi(mat_s)
mi_sw <- calculate_weighted_mi(phy_s,mat_s, spike_weights)


test <- apply(mat, 2, table )
test2 <- lapply(test, length)

test[[5]]
