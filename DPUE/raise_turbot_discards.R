### discards turbot for IBPTURBOT 2017
library(tidyverse)
library(stringr)
datapath <- "W:\\IMARES\\IJmuiden\\WOT\\WOT Discards demersaal\\3_data\\sasdata\\V201708_tarbot_discards\\"
output <- "W:\\IMARES\\Data\\ICES-WG\\IBPTURBOT\\2017\\DPUE\\data\\"
### got observer discards of turbot per lengthclass raised to triplevel from frisbe
### filter for TBB 70-99 and turbot discards
samples_obs <- 
  read_csv(paste(datapath, "nw_disc_trip_obs.csv", sep = ""), guess_max = 400000) %>%
  filter(str_detect(metier, "TBB"), str_detect(metier, "70"), DUTCH_NAME == "Tarbot") %>% 
  mutate(hpsegment2 = ifelse(POWER <= 221, "Euro", "Groot"),
         metier = if_else(metier == "TBB_DEF_70-99_0_0_G300hp","TBB_DEF_70-99", "TBB_DEF_70-99"),
         gear = if_else(str_detect(TOR_CODE, "e"), "pulse", "non-pulse")) ### only pulse and non-pulse because
                                                                          ### other distinctions are not clear in
write_csv(samples_obs, paste(output, "samples_obs.csv"))                  ### database

observer <- samples_obs %>%
  group_by(year, ship, week, quar, hpsegment2, gear, TRIPNR) %>%
  summarise(w_total = sum(w_total_trip), w_hour = sum(w_per_hour))

write_csv(observer, paste(output, "observer_aggr.csv"))

observer_trip <- observer %>% group_by(year) %>% 
  summarise(mean_w = mean(w_total))

write_csv(observer_trip, paste(output, "observer_year.csv"))

### got selfsampling discards of turbot per lengthclass raised to triplevel from frisbe 
samples_self <- 
  read_csv(paste(datapath, "nw_disc_trip_ss.csv", sep = ""), guess_max = 400000) %>%
  filter(str_detect(metier, "TBB"), str_detect(metier, "70"), DUTCH_NAME == "Tarbot") %>% 
  mutate(hpsegment2 = ifelse(POWER <= 221, "Euro", "Groot"),
         metier = if_else(metier == "TBB_DEF_70-99_0_0_G300hp","TBB_DEF_70-99", "TBB_DEF_70-99"),
         gear = if_else(str_detect(TOR_CODE, "e"), "pulse", "non-pulse"))

write_csv(samples_self, paste(output, "samples_self.csv"))

self <- samples_self %>%
  group_by(year, ship, week, quar, hpsegment2, gear, TRIPNR) %>%
  summarise(w_total = sum(w_total_trip), w_hour = sum(w_per_hour))

write_csv(self, paste(output, "self_aggr.csv"))

self_trip <- self %>% group_by(year) %>% 
  summarise(mean_w = mean(w_total))

write_csv(self_trip, paste(output, "self_year.csv"))

### got sample effort from visstat
effort_samples <-  
  read_csv(paste(datapath, "effort_trip_self.csv", sep = "")) %>%
  filter(metier == "TBB_DEF_70-99_0_0")

write_csv(effort_samples, paste(output, "effort_samples.csv"))

### join samples and effort
dpue_obs <-
  left_join(observer, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "observer")

dpue_obs_q <-
  left_join(observer, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, quar.x) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "observer")

dpue_obs_hp <-
  left_join(observer, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, hpsegment2.y) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "observer")

dpue_obs_hp_q <-
  left_join(observer, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, hpsegment2.y, quar.x) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "observer")

dpue_obs_hp_q_gear <-
  left_join(observer, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, hpsegment2.y, quar.x, gear) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "observer")

dpue_obs_q_gear <- 
  left_join(observer, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, quar.x, gear) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "observer")

dpue_obs_gear <- 
  left_join(observer, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, gear) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "observer")

dpue_self <-
  left_join(self, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "self")

dpue_self_q <-
  left_join(self, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, quar.x) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "self")

dpue_self_hp <-
  left_join(self, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, hpsegment2.y) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "self")

dpue_self_hp_q <-
  left_join(self, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, hpsegment2.y, quar.x) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "self")

dpue_self_hp_q_gear <-
  left_join(self, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, hpsegment2.y, quar.x, gear) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "self")

dpue_self_q_gear <- 
  left_join(self, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, quar.x, gear) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "self")

dpue_self_gear <- 
  left_join(self, effort_samples, by = "TRIPNR") %>% 
  mutate(wpue = w_total/hpeffort) %>% 
  group_by(year.x, gear) %>% 
  summarise(dpue = mean(wpue), n = n_distinct(TRIPNR)) %>% 
  mutate(prog = "self")

so <- full_join(dpue_self, dpue_obs)
write_csv(so, paste(output, "dpue.csv"))
soq <- full_join(dpue_self_q, dpue_obs_q) %>% 
  mutate(group = if_else(quar.x == "1",
                         if_else(prog == "self", "1_self", "1_observer"),
                         if_else(quar.x == "2",
                                 if_else(prog == "self", "2_self", "2_observer"),
                                  if_else(quar.x == "3",
                                          if_else(prog == "self", "3_self", "3_observer"),
                                          if_else(prog == "self", "4_self", "4_observer")))))
write_csv(soq, paste(output, "dpue_q.csv"))
                         
sog <- full_join(dpue_self_gear, dpue_obs_gear) %>% 
  mutate(group = if_else(gear == "pulse",
                          if_else(prog == "self", "pulse_self", "pulse_observer"),
                          if_else(prog == "self", "non-pulse_self", "non-pulse_observer")))
write_csv(sog, paste(output, "dpue_g.csv"))
# sohp <- full_join(dpue_self_hp, dpue_obs_hp) %>% 
#   mutate(grouper = if_else(hpsegment2.y == "Euro",
#                          if_else(prog == "self", "Euro_self", "Euro_observer"),
#                          if_else(prog == "self", "Groot_self", "Groot_observer")))
# sohpq <- full_join(dpue_self_hp_q, dpue_obs_hp_q) %>% 
#   mutate(group = if_else(hpsegment2.y == "Euro",
#                          if_else(prog == "self",
#                                  if_else(quar.x == 1, "Euro_self_1", 
#                                          if_else()),
#                          if_else(prog == "self", "Groot_self", "Groot_observer")))
# sogq <- full_join(dpue_self_q_gear, dpue_obs_q_gear)
# sohpqg <- full_join(dpue_self_hp_q_gear, dpue_obs_hp_q_gear)

### compare self vs observer
self_observer <- ggplot(so, aes(year.x, dpue, group = prog, linetype = prog)) +
  geom_line() +
  geom_text(aes(x = year.x, label = n)) +
  theme_bw()

quarter <- ggplot(soq, aes(year.x, dpue, group = group, linetype = prog)) + 
  geom_line(aes(colour = group)) +
  geom_text(aes(x = year.x, label = n)) +
  theme_bw()

gear <- ggplot(sog, aes(year.x, dpue, group = group, linetype = prog)) +
  geom_line(aes(colour = group)) +
  geom_text(aes(x = year.x, label = n)) +
  theme_bw()
