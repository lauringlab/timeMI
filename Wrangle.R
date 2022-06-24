#Wrangling alignment and metadata files
library(Biostrings)
library(tidyverse)
library(lubridate)

#Read in fasta files as AAStringSet
h3n2_pb2<- readAAStringSet("data/h3n2_pb2_1968_2021_aligned.fasta")
h3n2_pb1<- readAAStringSet("data/h3n2_pb1_1968_2021_aligned.fasta")
h3n2_pa<- readAAStringSet("data/h3n2_pa_1968_2021_aligned.fasta")
h3n2_np<- readAAStringSet("data/h3n2_np_1968_2021_aligned.fasta")
h3n2_ha<- readAAStringSet("data/h3n2_ha_1968_2021_aligned.fasta")

#At least one metadata is read in to maintain strain names
#Read in remaining if you want to do morestuff with dates
pb2_meta <- read_tsv("data/h3n2_pb2_1968_2021_f.tsv")
# pb1_meta <- read_tsv("data/h3n2_pb1_1968_2021_f.tsv")
# pa_meta <- read_tsv("data/h3n2_pa_1968_2021_f.tsv")
# np_meta <- read_tsv("data/h3n2_np_1968_2021_f.tsv")

#Function to turn dataframe back into AAStringSet
df_to_aastring <- function(df){
  aa <- AAStringSet(pull(df, seq))
  names(aa) <- df$rowname
  aa
}

#First, convert to dataframe
pb2_df <- as.data.frame(h3n2_pb2)
pb1_df <- as.data.frame(h3n2_pb1)
pa_df <- as.data.frame(h3n2_pa)
np_df <- as.data.frame(h3n2_np)
ha_df <- as.data.frame(h3n2_ha)

h3n2_ha_df <- full_join(rownames_to_column(pb2_df), rownames_to_column(pb1_df), by = "rowname") %>%
  full_join(rownames_to_column(pa_df), by = "rowname") %>%
  right_join(rownames_to_column(ha_df), by = "rowname") %>%
  unite("seq", x.x, x.y, x.x.x, x.y.y, sep = "") %>%
  left_join(pb2_meta %>% select(strain, date), by = c("rowname"="strain")) %>%
  #Remove weird czech sequence with a "J"
  filter(rowname != "A/Czech_Republic/216/2013///EPI_ISL_163601")

h3n2_ha_aa <- df_to_aastring(h3n2_ha_df)

#Save concatenated genes as single fasta
writeXStringSet(h3n2_ha_aa, "data/h3n2_pol_ha_aligned.fasta")

#Make seq_times dataframe
meta <- read_tsv("data/h3n2_pa_1968_2021_f.tsv")

seq_times <- as.data.frame(h3n2_ha_aa) %>%
  rownames_to_column() %>%
  left_join(meta, by = c("rowname" = "strain")) %>%
  mutate(date = year(date))

saveRDS(seq_times, "data/seq_times.rda")
