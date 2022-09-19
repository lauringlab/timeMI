library(tidyverse)
library(purrr)
library(data.table)
library(foreach)
library(doSNOW)
library(readxl)
library(lubridate)
library(Biostrings)

seq_times <- readRDS("data/seq_times_091122.rda")
seq_times_rndate <- seq_times %>%
  select(rowname, date)
h3n2_ha <- readAAStringSet("/Users/saraharcos/Desktop/H3N2_Processing/polha_aa.fasta") %>%
  as.matrix()

seq_times_noHA <- readRDS("data/seq_times_noHA_091122.rda")
h3n2_pol <- readAAStringSet("/Users/saraharcos/Desktop/H3N2_Processing/pol_aa.fasta") %>%
  as.matrix()

h3n2_mi <- calculate_mi(h3n2_ha)
saveRDS(h3n2_mi, "data/MI_results_091122/h3n2_mi.rds")

h3n2_mi_noHA <- calculate_mi(h3n2_pol)
saveRDS(h3n2_mi_noHA, "data/MI_results_091122/h3n2_mi_noHA.rds")


h3n2_meaned_mi_noHA <- calculate_meaned_mi(h3n2_pol, seq_times_noHA)
saveRDS(h3n2_meaned_mi_noHA, "data/MI_results_091122/h3n2_meaned_mi_noHA.rds")

h3n2_meaned_mi <- calculate_meaned_mi(h3n2_ha, seq_times)
saveRDS(h3n2_meaned_mi, "data/MI_results_091122/h3n2_meaned_mi.rds")


t <- paste(h3n2_ha[,1], h3n2_ha[,2])
names(t) <- rownames(h3n2_ha)
t3 <- t[seq_times_rndate[seq_times_rndate$date == 1968,]$rowname]

t1 <- h3n2_ha[seq_times_rndate[seq_times_rndate$date == 1968,]$rowname,]
t2 <- paste(t1[,1], t1[,2])
names(t2) <- rownames(t1)






run_mi <- function(phy, mat){
  print("calculating mi")
  mi <- calculate_mi(mat)
  print("calculating meaned mi")
  mi_m <- calculate_meaned_mi(mat, phy)
  return(list("mi" = mi, "mi_m" = mi_m))
}

calculate_mi <- function(mat){
  entropies <- vector(length = ncol(mat))
  for(i in 1:ncol(mat)){
    entropies[i] = get_ent_unweighted(i, mat)
  }
  
  non_zero_e <- which(entropies!=0)
  
  joint_cols <- t(combn(non_zero_e, 2)) %>%
    as.data.frame()
  
  # joint_ents <- vector(length = nrow(joint_cols))
  # for(i in 1:nrow(joint_cols)){
  #   joint_ents[i] = get_jointent_unweighted(joint_cols[i, "V1"], joint_cols[i, "V2"], mat)
  # }
  
  cl <- makeCluster(4, type="SOCK") # for 4 cores machine
  clusterExport(cl, c("get_jointent_unweighted", "calculate_mi"))
  registerDoSNOW (cl)
  
  pb <- txtProgressBar(max = nrow(joint_cols), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  # parallelization with vectorization
  joint_ents <- foreach(i = 1:nrow(joint_cols), .combine="c", .packages=c('tidyverse', 'data.table'),
                        .options.snow = opts) %dopar%
    {
      get_jointent_unweighted(joint_cols[i, "V1"], joint_cols[i, "V2"], mat)
    }
  
  close(pb)
  stopCluster(cl)
  
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
  
  return(sum(-t*log2(t)))
}

get_jointent_unweighted <- function(residue_a,residue_b, mat){
  jointcol <- as.matrix(paste(mat[,residue_a],mat[,residue_b],sep=""))
  
  t <- table(jointcol) %>%
    prop.table()
  
  return(sum(-t*log2(t)))
}

calculate_meaned_mi <- function(mat, df){
  ents <- c()
  for(i in 1:ncol(mat)){
    ents[i] = get_ent_meaned(i, df, mat)
    if(i%%100 == 0){
      print(paste(i, "out of", ncol(mat), sep  = " "))
    }
  }
  
  non_zero_e <- which(ents!=0)
  
  print("doing combn step")
  joint_cols <- t(combn(non_zero_e, 2)) %>%
    as.data.frame()

  print("starting clusters")
  cl <- makeCluster(6, type="SOCK") # for 4 cores machine
  clusterExport(cl, c("get_jointent_meaned", "calculate_meaned_mi"))
  registerDoSNOW (cl)

  pb <- txtProgressBar(max = nrow(joint_cols), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  # parallelization with vectorization
  joint_ents <- foreach(i = 1:nrow(joint_cols), .combine="c", .packages=c('tidyverse', 'data.table'),
                   .options.snow = opts) %dopar%
    {
      get_jointent_meaned(joint_cols[i, "V1"], joint_cols[i, "V2"], df, mat)
    }

  close(pb)
  stopCluster(cl)
  
  ents <- unlist(ents)
  return(format_mi(ents, joint_ents, joint_cols))
}

get_ent_meaned <- function(pos, df, mat){
  groups <- unique(df$date)
  
  t1f <- list()
  for(i in 1:length(groups)){
    # group_ids <- df %>%
    #   filter(date == groups[i]) %>%
    #   pull(rowname)
    # mat_i <- mat[group_ids,]
    mat_i <- mat[df[df$date == groups[i],]$rowname,]
    if(is.null(dim(mat_i))){
      t1f[[i]] = as.list(1)
    }
    else{
      t1f[[i]] <- as.list(prop.table(table(mat_i[,pos])))
    }
  }
  
  t1frbind <- rbindlist(t1f, fill = TRUE)
  t1frbind[is.na(t1frbind)] <- 0
  t1frw <- colMeans(t1frbind, na.rm = TRUE)
  
  return(sum(-log2(t1frw)*t1frw))
}

get_jointent_meaned <- function(a, b, df, mat){
  groups <- unique(df$date)
  
  # year_mat <- colmat[unlist(yeargeog_df[yeargeog_df$date_geog == year_geog,]$rownames),]
  
  t12f <- list()
  for(i in 1:length(groups)){
    # group_ids <- df %>%
    #   filter(date == groups[i]) %>%
    #   pull(rowname)
    # mat_i <- mat[group_ids,]
    mat_i <- mat[df[df$date == groups[i],]$rowname,]
    if(is.null(dim(mat_i))){
      t12f[[i]] = as.list(1)
    }
    else{
      t12f[[i]] <- as.list(prop.table(table(paste(mat_i[,a], mat_i[,b]))))
    }
  }
  
  t12frbind <- rbindlist(t12f, fill = TRUE)
  t12frbind[is.na(t12frbind)] <- 0
  t12frw <- colMeans(t12frbind, na.rm = TRUE)
  
  return(sum(-log2(t12frw)*t12frw))
}





