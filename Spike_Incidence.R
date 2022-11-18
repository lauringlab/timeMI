library(tidyverse)
library(readxl)
library(lubridate)


cases <- read_xlsx("/Users/saraharcos/Downloads/Cases and Deaths by County and by Date of Symptom Onset or by Date of Death2022-11-08.xlsx")
isolates <- read_csv("/Users/saraharcos/Downloads/wash_list.csv")

wash_cases <- cases %>%
  filter(COUNTY == "Washtenaw" &
           CASE_STATUS == "Confirmed" &
           Date > "2021-05-01" &
           Date < "2022-04-30") %>%
  mutate(Month_Yr = format_ISO8601(Date, precision = "ym")) %>%
  select(Month_Yr, Cases) %>%
  unique() %>%
  group_by(Month_Yr) %>%
  summarize(month_cases = sum(Cases))

wash_cases_isolates <- isolates %>%
  filter(coll_date > "2021-05-01" &
           coll_date <= "2022-04-30") %>%
  mutate(Month_Yr = format_ISO8601(coll_date, precision = "ym")) %>%
  count(Month_Yr) %>%
  full_join(wash_cases, by = c("Month_Yr")) %>%
  ungroup()


ggplot(wash_cases_isolates, aes(x = Month_Yr, y = n))+
  geom_col(fill = "grey") +
  geom_line(aes(y = month_cases/20,group = 1),
            size = 0.3) +
  geom_point(aes(y = month_cases/20),
             size = 0.5,
             alpha = 0.5) +
  xlab("Date") +
  ylab("Number of sequences") +
  labs(subtitle = "Distribution of complete H3N2 polymerase sequences in GISAID") +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 90))

ggsave("/Users/saraharcos/Desktop/MI Paper/Figures/pol_distribution.pdf",
       device = "pdf",
       height = 2.5, width = 3)








