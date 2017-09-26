rm(list = ls())
library(tidyverse)
datap <- "D:\\ICES_WG\\wg_IBPTur.27.4\\DPUE\\check\\"

dpue <- read_csv(paste(datap, "dpue_v2.csv", sep = ""))

ggplot(dpue) + geom_line(aes(year, dpue, group = prog, colour = prog))
