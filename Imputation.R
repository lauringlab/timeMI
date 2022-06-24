#Imputation
library(tidyverse)
library(data.table)
library(Biostrings)
library(purrr)
library(foreach)
library(doSNOW)


#First, look at distribution of gaps and X's
#Data
h3n2_ha <- readAAStringSet("data/h3n2_pol_ha_aligned.fasta")

h3n2_ha_mat <- h3n2_ha %>%
  as.matrix()

seq_times <- readRDS("data/seq_times.rda")

h3n2_x <- colSums(h3n2_ha_mat=="X")
h3n2_gap <- colSums(h3n2_ha_mat=="-")

x_plot <- h3n2_x %>%
  as.data.frame() %>%
  ggplot(aes(x = ./5038)) +
  geom_histogram(bins = 100) +
  labs(x = "Proportion of X residues", title = "Proportion of X residues per column")

saveRDS(x_plot, "x_distribution.rds")

gap_plot <- h3n2_gap %>%
  as.data.frame() %>%
  ggplot(aes(x = ./5038)) +
  geom_histogram(bins = 100) +
  labs(x = "Proportion of Gaps", title = "Proportion of Gaps per column")

saveRDS(gap_plot, "gap_distribution.rds")


h3n2_remove <- c(which(h3n2_x>length(h3n2_ha)/100), which(h3n2_gap>length(h3n2_ha)/100))

#DO NOT impute gaps. ONLY impute unknown (X) residues
h3n2_impute <- sum(h3n2_x > 0)

###############Removing above 1% since that is only 4 residues##################

################Now for the Imputation (by year and country)####################

######################## Weighting by Year and Continent #######################

meta <- read_tsv("/Users/saraharcos/Desktop/Lauring Lab/H3N2_Cooccurrence/data/gisaid_epiflu_isolates.txt") %>%
  select(Isolate_Id, Location) %>%
  mutate(Location = str_remove_all(Location, " ")) %>%
  mutate(Location =  str_extract(Location, "[A-z]*"))

seq_geog <- seq_times %>%
  left_join(meta, by = c("isolateID" =  "Isolate_Id")) %>%
  mutate(date_geog = paste(date, Location, sep = "_"))

#lookup tables for rowname:date and date:weight
lt_rowdategeog <- seq_geog %>%
  filter(!is.na(date)) %>%
  select(rowname, date_geog)

lt_weightsgeog <- lt_rowdategeog %>%
  count(date_geog) %>%
  mutate(inverse_weight = 1/n) %>%
  mutate(inverse_weight = inverse_weight/sum(inverse_weight)) %>%
  select(-n)


yeargeog_df <-lt_rowdategeog %>%
  nest(rownames = rowname) %>%
  full_join(lt_weightsgeog, by = "date_geog")

get_freqs_geog <- function(year_geog, impute = FALSE){
  print(paste("Running", year_geog, sep = " "))
  if (impute) {
    year_mat <- imputed_h3n2_mat[unlist(yeargeog_df[yeargeog_df$date_geog == year_geog,]$rownames),]
  }
  else{
    year_mat <- h3n2_ha_mat[unlist(yeargeog_df[yeargeog_df$date_geog == year_geog,]$rownames),]
  }
  
  year_weight <- lt_weightsgeog %>%
    filter(date_geog == year_geog) %>%
    pull(inverse_weight)
  
  if(is.null(dim(year_mat))){
    year_mat_freq <- lapply(year_mat, table)
  }
  else {
    year_mat_count <- apply(year_mat, 2, table, simplify = FALSE)
    year_mat_freq <- lapply(year_mat_count, prop.table)
  }
  return(lapply(year_mat_freq, function(x) as.list(x*year_weight)))
}

get_ent <- function(pos) {
  f <- colSums(pos, na.rm = TRUE)
  return(sum(-f*log(f)))
}

############################Finding Manual Gaps and Xs##############################
manual_x1 <- seq_geog %>%
  group_by(date_geog) %>%
  add_count(date_geog, name = "per_group") %>%
  summarize(q = sum(str_detect(x, "X")),
            per_group = per_group) %>%
  ungroup() %>%
  unique() %>%
  mutate(problem = q >= per_group)

manual_x <- seq_geog %>%
  filter(date_geog %in% c("1982_Oceania",
                          "1993_Oceania",
                          "1986_Europe",
                          "2013_Africa") &
           str_detect(x, "X"))

####################Manually replacing problematic X's##########################

h3n2_imputedx <- h3n2_ha
subseq(h3n2_imputedx["A/Victoria/186/1982//7000014_VW14/EPI_ISL_173279"], start = 2448, end = 2448) <- "R"
subseq(h3n2_imputedx["A/Sydney/2/1993/51285/93/9307027_VW1364/EPI_ISL_173280"], start = 2245, end = 2246) <- "LC"
subseq(h3n2_imputedx["A/Leningrad/360/1986/360//EPI_ISL_937"], start = 1693, end = 1693) <- "L"
subseq(h3n2_imputedx["A/Leningrad/360/1986/360//EPI_ISL_937"], start = 2389, end = 2389) <- "A"
subseq(h3n2_imputedx["A/Leningrad/360/1986/360//EPI_ISL_937"], start = 2437, end = 2437) <- "S"
subseq(h3n2_imputedx["A/Leningrad/360/1986/360//EPI_ISL_937"], start = 2471, end = 2471) <- "R"
subseq(h3n2_imputedx["A/Kenya/118/2013/A/KENYA/118/2013_ORIGINAL/2013761605/EPI_ISL_149681"], start = 461, end = 461) <- "V"
subseq(h3n2_imputedx["A/Kenya/118/2013/A/KENYA/118/2013_ORIGINAL/2013761605/EPI_ISL_149681"], start = 2445, end = 2445) <- "L"

# test <- r_x %>%
#   filter(date_geog == "2013_Africa") %>%
#   pull(x) %>%
#   str_split("", simplify = TRUE)
# 
# test2 <- seq_geog %>%
#   mutate(x = str_sub(x, start = 2472, end = 2472)) %>%
#   arrange(date_geog)

########################Programatic replacing X########################
h3n2_imp <- h3n2_imputedx %>%
  as.data.frame() %>%
  rownames_to_column()

impute <- function(year_geog){
  print(paste("Imputing: ", year_geog, sep = ""))
  observations <- h3n2_imp %>%
    filter(rowname %in% lt_rowdategeog[lt_rowdategeog$date_geog == year_geog,]$rowname)
  
  cols <- c()
  for (i in 1:str_length(observations[1,2])){
    col <- observations %>%
      mutate(x = str_sub(x, start = i, end = i)) %>%
      pull(x)
    
    if("X" %in% col) {
      col[which(col == "X")] <- names(sort(table(col),decreasing=TRUE)[1])
    }
    cols <- paste(cols, col, sep = "")
  }
  
  tibble(rowname = lt_rowdategeog[lt_rowdategeog$date_geog == year_geog,]$rowname,
         x = cols)
}

# t <- impute("2002_Oceania")

year_geogs <- unique(lt_weightsgeog$date_geog)

imputed_yeargeogs <- lapply(year_geogs, impute)
imputed_h3n2 <- rbindlist(imputed_yeargeogs)
imputed_h3n2_mat <- str_split_fixed(imputed_h3n2$x, pattern = "", n = Inf)
rownames(imputed_h3n2_mat) <- imputed_h3n2$rowname

saveRDS(imputed_h3n2_mat, "data/h3n2_pol_ha_aligned_imputed.rds")

##################Calculating Entropies w impute functions######################

ptm <- proc.time()
#Slow, 17 seconds
vgeog <- vector(mode = "list", length(year_geogs))
vgeog <- lapply(year_geogs, get_freqs_geog, impute = TRUE)

vlgeog <- vector(mode = "list", length(vgeog[[2]]))
vlgeog <- purrr::transpose(vgeog)

vmgeog <- vector(mode = "list", length(vlgeog))
vmgeog <- lapply(vlgeog, rbindlist, fill = TRUE)

vegeog <- vector(mode = "list", length(vmgeog))
vegeog <- lapply(vmgeog, get_ent)

#End time
proc.time() - ptm

#Add Joint entropies

non_zero_vegeog <- which(vegeog!=0)
non_zero_vegeog_removed <- non_zero_vegeog[!non_zero_vegeog %in% h3n2_remove]

joint_colsgeog <- t(combn(non_zero_vegeog_removed, 2)) %>%
  as.data.frame()

get_joint_freqs2geog <- function(year_geog, colmat){
  year_mat <- colmat[unlist(yeargeog_df[yeargeog_df$date_geog == year_geog,]$rownames),]
  
  year_mat_freq  <- year_mat %>% table() %>% prop.table()
  return(as.list(year_mat_freq*yeargeog_df[yeargeog_df$date_geog == year_geog,]$inverse_weight))
}

get_j_forgeog <- function(colA, colB, impute = FALSE){
  if (impute) {
    print("imputing")
    col_mat <- as.matrix(paste(imputed_h3n2_mat[,colA],imputed_h3n2_mat[,colB],sep=""))
    rownames(col_mat) <- rownames(imputed_h3n2_mat)
  }
  else{
    print("not imputing")
    col_mat <- as.matrix(paste(h3n2_ha_mat[,colA],h3n2_ha_mat[,colB],sep=""))
    rownames(col_mat) <- rownames(h3n2_ha_mat)
  }
  
  vjfor <- vector(mode = "list", length(year_geogs))
  for (i in 1:length(year_geogs)) {
    vjfor[[i]] = get_joint_freqs2geog(year_geogs[i], col_mat)
  }
  
  return(rbindlist(vjfor,fill = TRUE) %>% get_ent())
}

# parallel processing
cl <- makeCluster(4, type="SOCK") # for 4 cores machine
registerDoSNOW (cl)

pb <- txtProgressBar(max = nrow(joint_colsgeog), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
# parallelization with vectorization
tgeog <- foreach(i = 1:nrow(joint_colsgeog), .combine="c", .packages=c('tidyverse', 'data.table'),
                 .options.snow = opts) %dopar%
  {
    get_j_forgeog(joint_colsgeog[i, "V1"], joint_colsgeog[i, "V2"], impute = TRUE)
  }

joint_colsgeog$joint_entropy <- tgeog

close(pb)
stopCluster(cl) 

entropiesgeog <- unlist(vegeog)

v1_entropygeog <- numeric(nrow(joint_colsgeog))
v2_entropygeog <- numeric(nrow(joint_colsgeog))
for (i in 1:nrow(joint_colsgeog)) {
  v1_entropygeog[i] <- entropiesgeog[joint_colsgeog$V1[i]]
  v2_entropygeog[i] <- entropiesgeog[joint_colsgeog$V2[i]]
}

joint_colsgeog$v1_entropygeog <- v1_entropygeog
joint_colsgeog$v2_entropygeog <- v2_entropygeog

mi_hageog <- joint_colsgeog %>%
  mutate(mi = v1_entropygeog + v2_entropygeog - joint_entropy)

meanMIgeog <- mean(mi_hageog$mi)

mean_MIsgeog <- numeric(length(non_zero_vegeog_removed))

for (i in non_zero_vegeog_removed){
  mean_MIsgeog[i] <- mi_hageog[mi_hageog$V1 == i | mi_hageog$V2 == i,] %>% pull(mi) %>% mean()
}

v1_meanMIgeog <- numeric(nrow(joint_colsgeog))
v2_meanMIgeog <- numeric(nrow(joint_colsgeog))
for (i in 1:nrow(joint_colsgeog)) {
  v1_meanMIgeog[i] <- mean_MIsgeog[joint_colsgeog$V1[i]]
  v2_meanMIgeog[i] <- mean_MIsgeog[joint_colsgeog$V2[i]]
}

apc_hageog <- mi_hageog
apc_hageog$v1_meanMIgeog <- v1_meanMIgeog
apc_hageog$v2_meanMIgeog <- v2_meanMIgeog

saveRDS(apc_hageog, "/Users/saraharcos/Desktop/Lauring Lab/apc_hageog_imputed.rds")
saveRDS(apc_hageog, "data/apc_hageog_imputed.rds")
# apc_hageog <- readRDS("/Users/saraharcos/Desktop/Lauring Lab/apc_hageog_imputed.rds")
meanMIgeog <- mean(apc_hageog$mi)

MIp_hageog <- apc_hageog %>%
  mutate(apc = (v1_meanMIgeog*v2_meanMIgeog)/meanMIgeog,
         mip = mi - apc) %>%
  mutate(Group = paste("[", V1, ";", V2, "]", sep = ""))

saveRDS(MIp_hageog, "data/MIp_hageog_imputed.rds")
saveRDS(MIp_hageog, "/Users/saraharcos/Desktop/Lauring Lab/MIp_hageog_imputed.rds")


############################Run WITHOUT weighting###############################

get_ent_unweighted <- function(residue){
  t <- table(imputed_h3n2_mat[,residue]) %>%
    prop.table()
  
  return(sum(-t*log(t)))
}

entropies <- vector(length = ncol(imputed_h3n2_mat))
for(i in 1:ncol(imputed_h3n2_mat)){
  entropies[i] = get_ent_unweighted(i)
}

get_jointent_unweighted <- function(residue_a,residue_b){
  jointcol <- as.matrix(paste(imputed_h3n2_mat[,residue_a],imputed_h3n2_mat[,residue_b],sep=""))
  
  t <- table(jointcol) %>%
    prop.table()
  
  return(sum(-t*log(t)))
}

non_zero_e <- which(entropies!=0)
non_zero_e_removed <- non_zero_e[!non_zero_e %in% h3n2_remove]

joint_cols <- t(combn(non_zero_e_removed, 2)) %>%
  as.data.frame()

joint_ents <- vector(length = nrow(joint_cols))
for(i in 1:nrow(joint_cols)){
  joint_ents[i] = get_jointent_unweighted(joint_cols[i, "V1"], joint_cols[i, "V2"])
}

joint_cols$joint_entropy <- joint_ents

v1_entropy <- numeric(nrow(joint_cols))
v2_entropy <- numeric(nrow(joint_cols))
for (i in 1:nrow(joint_cols)) {
  v1_entropy[i] <- entropies[joint_cols$V1[i]]
  v2_entropy[i] <- entropies[joint_cols$V2[i]]
}

joint_cols$v1_entropy <- v1_entropy
joint_cols$v2_entropy <- v2_entropy

mi_ha <- joint_cols %>%
  mutate(mi = v1_entropy + v2_entropy - joint_entropy)

meanMI <- mean(mi_ha$mi)

mean_MIs <- numeric(length(non_zero_e_removed))

for (i in non_zero_e_removed){
  mean_MIs[i] <- mi_ha[mi_ha$V1 == i | mi_ha$V2 == i,] %>% pull(mi) %>% mean()
}

v1_meanMI <- numeric(nrow(joint_cols))
v2_meanMI <- numeric(nrow(joint_cols))
for (i in 1:nrow(joint_cols)) {
  v1_meanMI[i] <- mean_MIs[joint_cols$V1[i]]
  v2_meanMI[i] <- mean_MIs[joint_cols$V2[i]]
}

apc_ha <- mi_ha
apc_ha$v1_meanMI <- v1_meanMI
apc_ha$v2_meanMI <- v2_meanMI

saveRDS(apc_ha, "/Users/saraharcos/Desktop/Lauring Lab/apc_ha_imputed.rds")
saveRDS(apc_ha, "data/apc_ha_imputed.rds")

MIp_ha <- apc_ha %>%
  mutate(apc = (v1_meanMI*v2_meanMI)/meanMI,
         mip = mi - apc) %>%
  mutate(Group = paste("[", V1, ";", V2, "]", sep = ""))

saveRDS(MIp_ha, "data/MIp_ha_imputed.rds")
saveRDS(MIp_ha, "/Users/saraharcos/Desktop/Lauring Lab/MIp_ha_imputed.rds")







