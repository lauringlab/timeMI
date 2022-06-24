#Testing code with simulated data
#Simulated data generated using iqtree2 --alisim alignment_ph -t RANDOM{yh/7000} -m LG --length 2800

library(tidyverse)
library(purrr)
library(data.table)
library(foreach)
library(doSNOW)
library(readxl)
library(lubridate)

seq_times <- readRDS("/Users/saraharcos/Desktop/Lauring Lab/MI_Networks_App/data/seq_times.rda")
# seq_times <- readRDS("data/seq_times.rda")

sampling_weights <- seq_times %>%
  select(date) %>%
  filter(!is.na(date)) %>%
  count(date) %>%
  mutate(freq = n/sum(n)) %>%
  mutate(group = 1:47) %>%
  select(freq, group)


files <- list.files(path = "/Users/saraharcos/Downloads/iqtree-2.2.0-MacOSX/sims",
                    pattern = "phy",
                    full.names = TRUE)

# parallel processing
cl <- makeCluster(4, type="SOCK") # for 4 cores machine
registerDoSNOW (cl)

pb <- txtProgressBar(max = length(files), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
# parallelization with vectorization
foreach(i = 1:length(files), .combine="c", .packages=c('tidyverse', 'data.table'),
                 .options.snow = opts) %dopar%
  {
    run_sim(files[i])
  }

close(pb)
stopCluster(cl) 

run_sim <- function(file){
  filename <- str_extract(file, "[0-9]*.phy")
  
  phy <- read_phy(file)
  mat <- format_phy(phy)
  mi <- calculate_mi(mat)
  phy_s <- sample_phy(phy, sampling_weights)
  mat_s <- format_phy(phy_s)
  mi_s <- calculate_mi(mat_s)
  mi_sw <- calculate_weighted_mi(phy_s,mat_s)
  saveRDS(list("mi" = mi, "mi_s" = mi_s, "mi_sw" = mi_sw),
          paste("data/sim_mi/", filename, ".rds", sep = ""))
  return(list("mi" = mi, "mi_s" = mi_s, "mi_sw" = mi_sw))
}

phy <- read_phy(files[1])
mat <- format_phy(phy)
mi <- calculate_mi(mat)
phy_s <- sample_phy(phy, sampling_weights)
mat_s <- format_phy(phy_s)

test <- run_sim(files[1])

read_phy <- function(file){
  sim_data <- read_tsv(file) %>%
    mutate(`7000 100` = str_squish(`7000 100`)) %>%
    separate(col = `7000 100`, into = c("name", "seq"), sep = " ")
  return(sim_data)
}

format_phy <- function(phy){
  sim_mat <- c(str_split(phy$seq, "", simplify = TRUE))
  dim(sim_mat) <- c(7000,100)
  rownames(sim_mat) = phy$name
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
    seq = phy$seq[1:6956],
    group = sort(c(rep(1:47, 148)))
  ) %>%
    left_join(weights, by = "group") %>%
    slice_sample(n = 7000, replace = TRUE, weight_by = freq) %>%
    mutate(name = c(1:7000))
}

calculate_weighted_mi <- function(df, mat){
  inverse_weights <- df %>%
    select(name, group) %>%
    count(group) %>%
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


#############################Calculate Correlations#############################

sim_number <- numeric(1000)
sim_v_sample <- numeric(1000)
sim_v_wsample <- numeric(1000)
sample_v_wsample <- numeric(1000)

for (i in sim_number){
  #read rds
  #calculate correlation1 and save
  #calculate correlation2 and save
  #calculate correlation3 and save
}


sim <- readRDS("data/sim_mi/1.phy.RDS")

cor(sim$mi$mip, sim$mi_s$mip)
cor(sim$mi$mip, sim$mi_sw$mip)
cor(sim$mi_s$mip, sim$mi_sw$mip)

test <- sim$mi %>%
  left_join(sim$mi_s, by = "Group") %>%
  left_join(sim$mi_sw, by = "Group")

cor(test$mip.x, test$mip.y)

head(sim$mi_sw)

compare2 <- MIp %>%
  left_join(MIps, by = "Group") %>%
  left_join(MIp_hageog, by = "Group")

saveRDS(compare2, "simulation_compare.rds")

test <- readRDS("simulation_compare.rds")

ggplot(test, aes(x = mi.x, y = mi)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method='lm')


ggplot(test, aes(x = mip.x, y = mip.y)) +
  geom_point(alpha = 0.2) +
  # scale_x_log10() +
  # scale_y_log10() +
  geom_smooth(method='lm')

cor(test$mip.x, test$mip)
cor(test$mip.x, test$mip.y)


phy <- read_phy(files[1])
mat <- format_phy(phy)
mi <- calculate_mi(mat)
phy_s <- sample_phy(phy, sampling_weights)
mat_s <- format_phy(phy_s)
mi_s <- calculate_mi(mat_s)
mi_sw <- calculate_weighted_mi(phy_s,mat_s)





cor(mi$mip, mi_s$mip)
cor(mi$mip, mi_sw$mip)
  

  
  
  
  
  
  