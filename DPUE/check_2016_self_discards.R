rm(list = ls())
library(tidyverse)
datap <- "D:\\ICES_WG\\wg_IBPTur.27.4\\DPUE\\check\\"


year <- 2013 #change these each time for basefile
dyear <- 2009


effort <- read_csv(paste(datap, "\\", year, "\\effort_trip_self.csv", sep = "")) %>% mutate(kWdays = hpeffort * 0.745699872) %>% 
  filter(year == dyear)


self <- read_csv(paste(datap, "\\", year, "\\nw_disc_trip_ss.csv", sep = ""), guess_max = 40000) %>%
  filter(SCIENTIFIC_NAME == "Scophthalmus maximus", prog == "self sampling", year == dyear) %>%
  mutate(tripnr = TRIPNR + 0) %>% 
  group_by(tripnr, metier) %>% 
  summarise(w_total = sum(w_total_trip)) %>% 
  ungroup()

selfs <- left_join(self, effort, by = "tripnr") %>% filter(metier == "TBB_DEF_70-99_0_0") %>% mutate(dpue = w_total / kWdays)

dpue <- summarise(selfs, dpue = mean(dpue), n = n_distinct(tripnr)) %>% mutate(prog = "self")


ob <- read_csv(paste(datap, "\\", year, "\\nw_disc_trip_obs.csv", sep = ""), guess_max = 40000) %>%
  filter(SCIENTIFIC_NAME == "Scophthalmus maximus", prog == "self sampling", year == dyear) %>%
  mutate(tripnr = TRIPNR + 0) %>% 
  group_by(tripnr, metier) %>% 
  summarise(w_total = sum(w_total_trip)) %>% 
  ungroup()

obs <- left_join(ob, effort, by = "tripnr") %>% filter(metier == "TBB_DEF_70-99_0_0") %>% mutate(dpue = w_total / kWdays)
 
dpue_ob <- summarise(obs, dpue = mean(dpue), n = n_distinct(tripnr)) %>% mutate(prog = "obs")


dpue <- full_join(dpue, dpue_ob)

write_csv(dpue, paste(datap, "\\", dyear, "dpue.csv", sep = ""))
