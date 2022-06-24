#Scratch
library(Biostrings)
library(tidyr)
library(stringr)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(plotly)

h3n2_ha<- readAAStringSet("/Users/saraharcos/Desktop/Lauring Lab/H3N2_Cooccurrence/data/h3n2_pol_ha_aligned.fasta")
seq_times_ex <- readRDS("/Users/saraharcos/Desktop/Lauring Lab/MI_Networks_App/data/seq_times.rda") %>%
  select(rowname, date) %>%
  filter(!is.na(date))

#If I store the indices and pairs that I care about, this takes less memory than re-storing all the sequences

#This function inputs an MSA and outputs a formatted MSA, with each character it's own column
format_MSA <- function(msa){
  msa %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column()
}

t <- format_MSA(h3n2_ha)

#This function inputs a formatted MSA and outputs a vector of non-conserved *indices*
remove_conserved <- function(formatted_msa){
  unique_count <- formatted_msa %>%
    summarize_all(n_distinct) %>%
    select_if(function(col) col > 1)
  
  nonconserved <- colnames(unique_count) %>%
    str_remove("_1")
  
  nox <- t %>%
    select(-rowname) %>%
    summarize_all(function(c) c == "X") %>%
    summarize_all(sum) %>%
    select_if(function(col) col == 0)
  
  noxnames <- colnames(nox)
  
  nod<- t %>%
    select(-rowname) %>%
    summarize_all(function(c) c == "-") %>%
    summarize_all(sum) %>%
    select_if(function(col) col == 0)
  
  nodnames <- colnames(nod)
  
  intersect(noxnames, nonconserved) %>% intersect(nodnames)
}

m <- remove_conserved(t)

#This function takes a dataframe of observations and dates and adds a weight column based on the year
add_weights <- function(seq_times){
  n_years <- seq_times %>%
    count(date) %>%
    mutate(inverse_weight = 1/n) %>%
    mutate(inverse_weight = inverse_weight/sum(inverse_weight)) %>%
    select(-n)
  
  seq_times %>%
    left_join(n_years, by = "date")
}

weights <- add_weights(seq_times_ex)

get_col <- function(ind, formatted_msa) {
  formatted_msa %>%
    select("index" = all_of(ind), rowname)
}

a <- get_col("V4", t)



get_joint_col <- function(index1, index2, formatted_msa){
  formatted_msa %>%
    select(index1 = index1, index2 = index2, rowname) %>%
    mutate(index = paste0(index1, index2)) %>%
    select(index, rowname)
}

j <- get_joint_col("V3", "V4", t)

get_entropy <- function(col, weights){
  col %>%
    left_join(weights, by = "rowname") %>%
    select(-rowname) %>%
    filter(!is.na(date)) %>%
    group_by(date,inverse_weight, index) %>%
    tally() %>%
    mutate(freq = n/sum(n)) %>%
    mutate(weighted_freq = freq*inverse_weight) %>%
    ungroup() %>%
    group_by(index) %>%
    summarize(weighted_mean = sum(weighted_freq)) %>%
    mutate(freqlog = weighted_mean*log(weighted_mean)) %>%
    summarize(ent  = -sum(freqlog))
}

test <- get_entropy(a, weights)

get_mi <- function(index1, index2, weights, formatted_msa){
  index1_ent <- get_entropy(get_col(index1, formatted_msa), weights)
  index2_ent <- get_entropy(get_col(index2, formatted_msa), weights)
  joint_ent <- get_entropy(get_joint_col(index1, index2, formatted_msa), weights)
  
  return((index1_ent + index2_ent - joint_ent)[1,1])
}


get_mi("V3", "V4", weights,t)


get_mi_v <- Vectorize(get_mi, vectorize.args = c("index1", "index2"))

test <- t(combn(m, 2)) %>%
  as.data.frame()

test10 <- test[1:200,]
  

system.time(
  test_mi_100 <- test10 %>%
  mutate(mi = get_mi_v(V1, V2, weights, t))
  )

#saveRDS(test_mi, "/Users/saraharcos/Desktop/Lauring Lab/weighted_mi.RDS")

test_mi <- readRDS("/Users/saraharcos/Desktop/Lauring Lab/weighted_mi.RDS")

test_mi_1000 <- test_mi %>%
  mutate(mi = unlist(unname(mi))) %>%
  top_n(5000, mi) %>%
  mutate(Group = paste(V1, V2, sep = ";")) %>%
  mutate(Group = str_remove_all(Group, pattern = "V"))

apc_iav <- read_tsv("/Users/saraharcos/Desktop/Lauring Lab/H3N2_Cooccurrence/data/mica_results_HA.tsv") %>%
  mutate(Group = str_extract(string = Group, pattern = "[0-9]*\\;[0-9]*")) %>%
  separate(Group, into = c("a", "b"), sep = ";", remove = FALSE) %>%
  mutate(a = as.numeric(a),
         b = as.numeric(b))


apc_iav_1000 <- apc_iav %>%
  top_n(5000, MI) %>%
  select(Group, MI)


mica_compare <- test_mi_1000 %>%
  left_join(apc_iav_1000, by = "Group")

pb2 <- readRDS("/Users/saraharcos/Desktop/Lauring Lab/H3N2_Cooccurrence/data/pb2_length.rds")
pb1 <- readRDS("/Users/saraharcos/Desktop/Lauring Lab/H3N2_Cooccurrence/data/pb1_length.rds")
pa <- readRDS("/Users/saraharcos/Desktop/Lauring Lab/H3N2_Cooccurrence/data/pa_length.rds")

adjust_index <- function(index){
  
  domain = "none"
  
  if(index <= pb2){
    domain = "PB2"
  }
  else if(index <= pb2+pb1){
    domain = "PB1"
  }
  else if(index <= pb2+pb1+pa){
    domain = "PA"
  }
  else{
    domain = "HA"
  }
  
  new_index = 0
  
  if(domain == "PB2"){
    new_index = index
  }
  else if(domain == "PB1"){
    new_index = as.numeric(index - pb2)
  }
  else if(domain == "PA"){
    new_index = as.numeric(index - (pb2+pb1))
  }
  else{
    new_index = as.numeric(index - (pb2+pb1+pa))
  }
  
  return(paste(domain, new_index, sep = "_"))
}

adjust_index(3000)

adjust_index_v = Vectorize(adjust_index)

mica_compare <- test_mi_1000 %>%
  left_join(apc_iav_1000, by = "Group") %>%
  separate(Group, into = c("a", "b")) %>%
  mutate(a = adjust_index_v(as.numeric(a)),
         b = adjust_index_v(as.numeric(b))) %>%
  mutate(Group = paste(a,b, sep = ";"))


plog<- ggplot(mica_compare, aes(x = log(mi), y = log(MI), group = Group)) +
  geom_point(color = "steelblue", alpha = 0.5) +
  labs(x = "log10(weighted MI)", y = "log10(MI)",
       title = "Weighted MI versus MI, log transformed")

p <- ggplot(mica_compare, aes(x = mi, y = MI, group = Group)) +
  geom_point(color = "steelblue", alpha = 0.5) +
  labs(x = "weighted MI", y = "MI",
       title = "Weighted MI versus MI")

saveRDS(p, "/Users/saraharcos/Desktop/Lauring Lab/Lab_Meetings/weightedMIplot")
saveRDS(plog, "/Users/saraharcos/Desktop/Lauring Lab/Lab_Meetings/weightedMIlogplot")

ggplotly(plog)

ggplotly(p)

ta <- a %>% left_join(weights, by = "rowname") %>%
  select(-rowname) %>%
  filter(!is.na(date)) %>%
  group_by(date,inverse_weight, index) %>%
  tally() %>%
  mutate(freq = n/sum(n)) %>%
  mutate(weighted_freq = freq*inverse_weight) %>%
  ungroup() %>%
  group_by(index) %>%
  summarize(weighted_mean = sum(weighted_freq)) %>%
  mutate(freqlog = weighted_mean*log(weighted_mean)) %>%
  summarize(ent  = -sum(freqlog))

