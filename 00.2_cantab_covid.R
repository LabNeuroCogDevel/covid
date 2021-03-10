#!/usr/bin/Rscript

# combine cantab MOT/SOC summaries and covid questioneer
# 
# 20210310 WF

library(dplyr)
library(tidyr)
library(lubridate)


# TODO: get covid score instead of all of covid survey?
covid <- read.csv('data/qualtrics_adult-kid.csv') %>%
    rename(ID=External.Data.Reference) %>%
    mutate(vdate=gsub(' .*', '', Start.Date) %>% ymd()) %>%
    filter(!is.na(ID))

cantab <- read.csv('./data/cantab_pet7T_motsoc.csv') %>%
    mutate(vdate=ymd(vdate))

covid_cantab <- LNCDR::date_match(covid, cantab, 'ID', 'vdate', maxdatediff=Inf)

nrow(covid)        # 117
nrow(covid_cantab) # 117

write.csv(covid_cantab, 'data/covid_cantab.csv', row.names=F)
