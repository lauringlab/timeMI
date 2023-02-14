library(tidyverse)
library(readxl)
library(Biostrings)
library(lubridate)

# Data accessed from GISAID on August 30th, 2022
# GISAID Options:
#   Type: A
#   H   : 3
#   N   : 2
#   Host: Human
#   Loc : All
#   Required Segments: PB2, PB1, PA, HA
#   Only Complete: Yes

#Result: 7250 sequences

################################################################################

#Read in metadata from GISAID
mt <- read_xls("../H3N2_Processing/gisaid_20220830.xls")

# Identify sequences that are not the correct length
pb2_aa <- readAAStringSet("../H3N2_Processing/pb2_aa.fasta")
pb2_nfl <- names(pb2_aa[which(width(pb2_aa) != 759)])

pb1_aa <- readAAStringSet("../H3N2_Processing/pb1_aa.fasta")
pb1_nfl <- names(pb1_aa[which(width(pb1_aa) != 757)])

pa_aa <- readAAStringSet("../H3N2_Processing/pa_aa.fasta")
pa_nfl <- names(pa_aa[which(width(pa_aa) != 716)])

ha_aa <- readAAStringSet("../H3N2_Processing/ha_aa.fasta")
ha_nfl <- names(ha_aa[which(width(ha_aa) != 566)])

#Perform filtering
mt_filtered <- mt %>%
  #Remove egg passage
  filter(!str_detect(Passage_History, "egg|Egg|E[0-9]|AM|Am|E") |
           is.na(Passage_History)) %>%
  #Remove duplicate entries
  filter(!str_detect(`PB1 Segment_Id`, ",") &
           !str_detect(`HA Segment_Id`, ",")) %>%
  #Remove incorrect length
  filter(!Isolate_Id %in% c(pb2_nfl, pb1_nfl, pa_nfl, ha_nfl))

# Result: 6936 sequences

#Filter AA fasta files and save for alignment
writeXStringSet(pb2_aa[mt_filtered$Isolate_Id],
                "../H3N2_Processing/pb2_aa_filt.fasta")
writeXStringSet(pb1_aa[mt_filtered$Isolate_Id],
                "../H3N2_Processing/pb1_aa_filt.fasta")
writeXStringSet(pa_aa[mt_filtered$Isolate_Id],
                "../H3N2_Processing/pa_aa_filt.fasta")
writeXStringSet(ha_aa[mt_filtered$Isolate_Id],
                "../H3N2_Processing/ha_aa_filt.fasta")

#Concatenate aligned genes
#Aligned using MAFFT v7.490 (2021/Oct/30)
#Order of genes: PB2 - PB1 - PA - HA

pb2_aa_al <- readAAStringSet("../H3N2_Processing/pb2_aa_al.fasta")
pb1_aa_al <- readAAStringSet("../H3N2_Processing/pb1_aa_al.fasta")
pa_aa_al <- readAAStringSet("../H3N2_Processing/pa_aa_al.fasta")
ha_aa_al <- readAAStringSet("../H3N2_Processing/ha_aa_al.fasta")

#Polymerase only
pol_aa <- AAStringSet()
for(i in 1:length(pb2_aa_al)){
  pol_aa[i] <- xscat(pb2_aa_al[i], pb1_aa_al[i], pa_aa_al[i])
}
names(pol_aa) <- names(pb2_aa_al)

writeXStringSet(pol_aa, "../H3N2_Processing/pol_aa.fasta")
pol_aa <- readAAStringSet("../H3N2_Processing/pol_aa.fasta")

#Polymerase and HA
polha_aa <- AAStringSet()
for(i in 1:length(pb2_aa_al)){
  polha_aa[i] <- xscat(pb2_aa_al[i], pb1_aa_al[i], pa_aa_al[i], ha_aa_al[i])
}
names(polha_aa) <- names(pb2_aa_al)

writeXStringSet(polha_aa, "../H3N2_Processing/polha_aa.fasta")
polha_aa <- readAAStringSet("../H3N2_Processing/polha_aa.fasta")



########################## Make seq_times dataframes ##########################

seq_times <- as.data.frame(polha_aa) %>%
  rownames_to_column() %>%
  left_join(mt_filtered, by = c("rowname" = "Isolate_Id")) %>%
  mutate(date = year(parse_date_time(Collection_Date, c("ymd", "ym", "y"))))  %>%
  filter(!is.na(date))

saveRDS(seq_times, "data/seq_times_091122.rda")
saveRDS(seq_times, "/Users/sarcos/Desktop/Lauring Lab/MI_Networks_App/data/seq_times_091122.rda")

seq_times_noHA <- as.data.frame(pol_aa) %>%
  rownames_to_column() %>%
  left_join(mt_filtered, by = c("rowname" = "Isolate_Id")) %>%
  mutate(date = year(parse_date_time(Collection_Date, c("ymd", "ym", "y"))))  %>%
  filter(!is.na(date))

saveRDS(seq_times_noHA, "data/seq_times_noHA_091122.rda")


######################## Make mi_network_info #################################
mi_network_info <- tibble(
  index = 1:2798,
  gene = c(rep("PB2", 759), rep("PB1", 757), rep("PA", 716), rep("HA", 566)),
  id = c(1:759, 1:757, 1:716, -16:-1, 1:550)) %>%
  mutate(id = paste(gene, id, sep = "_"),
         color = case_when(
           str_detect(id, "PA") ~ "darkgrey",
           str_detect(id, "PB1") ~ "#08519c",
           str_detect(id, "PB2") ~ "#e41a1c",
           TRUE ~ "#fd8d3c"
         )) %>%
  select(-gene)

#antigenic site reference: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0092153#pone.0092153.s004
siteA <- c(122, 124, 126, 130:133,
           135, 137, 138, 140, 142:146,
           150, 152, 168)
siteB <- c(128, 129, 155:160, 163:165,
           186:190, 192:194, 196:198)
siteC <- c(44:48, 50, 51, 53, 54, 273,
           275, 276, 278:280, 294, 297,
           299, 300, 304, 305, 307:312)
siteD <- c(96, 102, 103, 117, 121, 167,
           170:177, 179, 182, 201, 203,
           207:209, 212:219, 226:230, 238,
           240, 242, 244, 246:248)
siteE <- c(57, 59, 62, 63, 67, 75, 78,
           80:83, 86:88, 91, 92, 94,
           109, 260:262, 265)

mi_network_info_ha <- tibble(
  index = 2233:2798,
  gene = c(rep("HA", 566)),
  id = c(-16:-1, 1:550)) %>%
  mutate(color = case_when(
           id %in% siteA ~ "lightpink",
           id %in% siteB ~ "darkblue",
           id %in% siteC ~ "lightgreen",
           id %in% siteD ~ "orange",
           id %in% siteE ~ "lightblue",
           TRUE ~ "lightgrey"
         ),
         id = paste(gene, id, sep = "_")) %>%
  select(-gene)

saveRDS(mi_network_info, "data/mi_network_info.rds")
saveRDS(mi_network_info, "/Users/sarcos/Desktop/Lauring Lab/MI_Networks_App/data/mi_network_info.rds")

saveRDS(mi_network_info_ha, "data/mi_network_info_ha.rds")
saveRDS(mi_network_info_ha, "/Users/sarcos/Desktop/Lauring Lab/MI_Networks_App/data/mi_network_info_ha.rds")


