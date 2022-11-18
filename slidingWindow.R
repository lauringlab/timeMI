#Run MI and Entropy on sliding 5-year windows across dataset, no HA***
#Don't need to calculate meaned or weighted MI

library(tidyverse)
library(purrr)
library(data.table)
library(foreach)
library(doSNOW)
library(readxl)
library(lubridate)

seq_times <- readRDS("data/seq_times_091122.rda") %>%
  filter(!is.na(date))
imputed_h3n2_mat<- readRDS("data/h3n2_pol_ha_aligned_imputed.rds")

start_years <- seq(min(seq_times$date), (max(seq_times$date) - 4))

window_results <- vector(mode = "list")
saveRDS(window_results, "data/slidingWindows.rds")
for (i in start_years){
  window_results[[i]] <- get_window(seq_times, imputed_h3n2_mat, i) %>%
    calculate_mi()
  print(paste("Done with window at", i))
}



get_window <- function(df, mat, start_year){
  window_df <- df %>%
    filter(date >= start_year & date < start_year+5)
  return(mat[window_df$rowname,])
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
  
  cl <- makeCluster(8, type="SOCK") # for 4 cores machine
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
  
  return(sum(-t*log(t)))
}

get_jointent_unweighted <- function(residue_a,residue_b, mat){
  jointcol <- as.matrix(paste(mat[,residue_a],mat[,residue_b],sep=""))
  
  t <- table(jointcol) %>%
    prop.table()
  
  return(sum(-t*log(t)))
}




windows <- rbindlist(window_results, idcol = TRUE)

windows_top <- windows %>%
  slice_max(n = 10000, order_by = mi)

ggplot(windows_top, aes(x = .id, y = mi, group = Group)) +
  geom_line(alpha = 0.3) +
  theme_minimal()

windows_var <- windows %>%
  group_by(Group) %>%
  mutate(mi_var = var(mi)) %>%
  select(Group, mi_var) %>%
  unique()

windows_var_top <- windows_top %>%
  group_by(Group) %>%
  mutate(mi_var = var(mi)) %>%
  select(Group, mi_var) %>%
  unique()

windows_med <- windows %>%
  group_by(Group) %>%
  mutate(mi_med = median(mi)) %>%
  select(Group, mi_med)

windows_med_top <- windows_med %>%
  unique()

test <- windows_med_top %>%
  ungroup() %>%
  slice_max(n=100, order_by = mi_med, with_ties = FALSE)

windows_filt <- windows %>%
  filter(Group %in% test$Group)

ggplot(windows_filt, aes(x = .id, y = mi, group = Group)) +
  geom_line(alpha = 0.3) +
  theme_minimal()

window_examples <- windows_filt %>%
  filter(Group %in% c("[1229;1850]", "[939;1346]", "[522;1846]", "[221;293]")) %>%
  select(.id, Group, mi, mip, v1_entropy, v2_entropy) %>%
  full_join(tibble(.id = start_years), by = ".id") %>%
  complete(.id, Group, fill = list(mi = 0, mip = -Inf, v1_entropy = 0, v2_entropy = 0)) %>%
  filter(!is.na(Group))

ggplot(window_examples, aes(x = .id, y = v1_entropy, group = Group, color = Group)) +
  geom_line() +
  theme_minimal()

get_position_freq <- function(position, seqs){
  #takes a position and gets the frequency of residues
  temp <- mutate(seqs, `Amino Acid` = str_sub(x, start = position, end = position)) %>%
    count(`Amino Acid`, date, .drop = FALSE) %>%
    group_by(date) %>%
    mutate(freq = n/sum(n)) %>%
    ungroup() %>%
    complete(date, `Amino Acid`, fill = list(n = 0, freq = 0))
  temp
}

get_subgraph_freqs <- function(elements){
  names(elements) <- elements
  m <- lapply(elements, get_position_freq, seq_times) %>%
    bind_rows(.id = "position") %>%
    mutate(position = get_index_v(position))
  return(m)
}

pal <- c("#a6cee3",
         "#1f78b4",
         "#b2df8a",
         "#33a02c",
         "#fb9a99",
         "#e31a1c",
         "#fdbf6f",
         "#ff7f00",
         "#cab2d6",
         "#6a3d9a",
         "#ffff99",
         "#b15928",
         "#bebebe",
         "#d2b48c",
         "#fccde5",
         "#ccebc5",
         "#fddbc7", 
         "#d4af37",
         "white",
         "darkblue",
         "turquoise")

pal2 <- c("A" = "#8CFF8C",
          "G" = "#FFFFFF",
          "L" = "#455E45",
          "S" = "#FF7042",
          "V" = "#FF8CFF",
          "T" = "#B84C00",
          "K" = "#4747B8",
          "D" = "#A00042",
          "I" = "#004C00",
          "N" = "#FF7C70",
          "E" = "#660000",
          "P" = "#525252",
          "R" = "#00007C",
          "F" = "#543C42",
          "Q" = "#FF4C4C",
          "Y" = "#8C704C",
          "H" = "#7070FF", 
          "C" = "#FFFF70",
          "M" = "#B8A042",
          "W" = "#4F4600",
          "X" = "#B8B8B8")

mi_iav_info <- readRDS("data/mi_network_info.rds")
get_index <- function(pos){
  return(mi_iav_info[mi_iav_info$index == pos,]$id %>%
           str_replace("_", " "))
}
get_index_v <- Vectorize(get_index)

plot_example_freqs <- function(v1, v2){
  ggplot(get_subgraph_freqs(c(v1, v2)), aes(x = date, y = freq, fill = `Amino Acid`,
                                            color = `Amino Acid`, group = `Amino Acid`)) +
    geom_area(alpha = 0.6) +
    facet_wrap(~position, ncol = 1) +
    scale_fill_manual(values = pal2) +
    scale_color_manual(values = pal2) +
    scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
    theme_minimal() +
    labs(title = "Amino Acid Frequencies",
         subtitle = paste(get_index(v1), "and", get_index(v2)),
         x = "Date", y = "Frequency") +
    theme(text = element_text(size = 12))
}

plot_example_mi <- function(v1, v2){
  
  window_example <- windows_filt_max %>%
    filter(V1 == v1 & V2 == v2) %>%
    select(.id, Group, mi, mip, v1_entropy, v2_entropy) %>%
    full_join(tibble(.id = start_years), by = ".id") %>%
    complete(.id, Group, fill = list(mi = 0, mip = -Inf, v1_entropy = 0, v2_entropy = 0)) %>%
    filter(!is.na(Group))
  
  ggplot(window_example, aes(x = .id, y = mi)) +
    geom_line() +
    scale_y_continuous(breaks = c(0, 0.3, 0.6)) +
    theme_minimal() +
    labs(title = "Mutual Information (5-year sliding window)",
         subtitle = paste(get_index(v1), "and", get_index(v2)),
         x = "Date", y = "Mutual Information") +
    theme(text = element_text(size = 12))
}

plot_example_freqs(939, 1346)
ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/freq_939_1346.pdf",
       width = 4, height = 3)
plot_example_mi(939, 1346)
ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/mi_939_1346.pdf",
       width = 4, height = 3)

plot_example_freqs(1229, 1850)
ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/freq_1220_1850.pdf",
       width = 4, height = 3)
plot_example_mi(1229, 1850)
ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/mi_1229_1850.pdf",
       width = 2.55, height = 1.825)

plot_example_freqs(227, 939)
ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/freq_227_939.pdf",
       width = 4, height = 3)
plot_example_mi(227, 939)
ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/mi_227_939.pdf",
       width = 2.55, height = 1.825)

plot_example_freqs_multi <- function(nodes){
  ggplot(get_subgraph_freqs(nodes), aes(x = date, y = freq, fill = `Amino Acid`,
                                            color = `Amino Acid`, group = `Amino Acid`)) +
    geom_area(alpha = 0.6) +
    facet_wrap(~position, ncol = 1) +
    scale_fill_manual(values = pal2) +
    scale_color_manual(values = pal2) +
    scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
    theme_minimal() +
    labs(title = "Amino Acid Frequencies",
         x = "Date", y = "Frequency") +
    theme(text = element_text(size = 12))
}

pb2 <- 759
pb1 <- 757
pa <- 716

plot_example_freqs_multi(c(227, 338, 194, 569))

ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/freq_subgraph1.pdf",
       width = 4, height = 3)

plot_example_freqs_multi(c(559, 697,
                           (573+pb2+pb1),
                           (343+pb2+pb1),
                           (557+pb2+pb1),
                           (312+pb2+pb1)))

ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/freq_subgraph2.pdf",
       width = 4, height = 4)


plot_example_freqs_multi(c(569, 338, 194,
                           (11+pb2+pb1+pa)))

ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/freqHA_subgraph5.pdf",
       width = 4, height = 3)

plot_example_freqs_multi(c(697, (343+pb2+pb1),
                           (573+pb2+pb1),
                           (312+pb2+pb1),
                           (73+pb2+pb1+pa)))

ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/freqHA_subgraph9.pdf",
       width = 4, height = 3.5)

plot_example_freqs_multi(c(44, 199, 591, 627, 645,
                           (321+pb2+pb1)))

ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/freq_627.pdf",
       width = 4, height = 4)

plot_example_freqs_multi(c((619+pb2), (709+pb2)))

ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/freq_subgraph4.pdf",
       width = 4, height = 2)

plot_example_freqs_multi(c((52+pb2), (576+pb2)))

ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/freq_subgraph8.pdf",
       width = 4, height = 2)

plot_example_freqs_multi(c(107,
                           (469+pb2),
                           (350+pb2+pb1)))

ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/freq_subgraph10.pdf",
       width = 4, height = 2.5)





windows_max <- windows %>%
  group_by(Group) %>%
  mutate(mi_med = max(mi)) %>%
  select(Group, mi_med) %>%
  unique() %>%
  ungroup() %>%
  slice_max(n=300, order_by = mi_med, with_ties = FALSE)

windows_filt_max <- windows %>%
  filter(Group %in% windows_max$Group)





ggplot(windows_var, aes(x = mi_var)) +
  geom_density() +
  scale_x_log10()

ggplot(windows_var_top, aes(x = mi_var)) +
  geom_density() +
  scale_x_log10()

ggplot(windows, aes(x = mi)) +
  geom_density(fill = "darkgreen", alpha = 0.3,
               color = "darkgreen") +
  scale_x_log10()

