library(tidyverse)
library(Biostrings)
library(lubridate)

sc2_full <- readDNAStringSet("data/gisaid_hcov-19_2022_06_16_22.fasta")
spike <- sc2_full[width(sc2_full) == 29782] %>% subseq(start = 21509, end = 25330)

spikec <- spike[names(spike) != "hCoV-19/USA/MI-UM-10036994331/2020|EPI_ISL_913517|2020-11-17" &
                  names(spike) != "hCoV-19/USA/MI-UM-S768053/2021|EPI_ISL_1585913|2021-03-27"]

spiken <- vcountPattern("N", spikec)
spikek <- vcountPattern("K", spikec)
spikex <- spikec[spiken == 0 & spikek == 0]
translate(spikex)

writeXStringSet(spikex, "data/spikent.fasta")

spikent_aligned <- readDNAStringSet("data/spikent_aligned.fasta")

spikep <- translate(spikent_aligned)

spike_df <- as.data.frame(spikep) %>%
  rownames_to_column() %>%
  mutate(coldate = str_sub(rowname, start = -10),
         coldate = as_date(coldate),
         date = format_ISO8601(coldate, precision = "ym"))


ggplot(spike_df, aes(x = colmonth)) +
  geom_bar()

###Need to make: fasta, seqtimes

saveRDS(spike_df, "data/spike_times.rds")




