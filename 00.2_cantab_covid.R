#!/usr/bin/Rscript

# combine cantab MOT+SOC summaries and covid questioneer

# uses                          | from
# ------------------------------|-------------------
#  data/cantab_pet7T_motsoc.csv | 00.1_cantab_merge.R 
#  data/qualtrics_adult-kid.csv | qualtrics_covid.R
#
# NB.
#  likely to have many CANTAB visits, but only one covid survey per subject
#  (117 uniq IDs as of 20210310)
#  the covid survey and CANTAB task could be any duration apart
#  using maxdatediff=Inf to allow matches that are years apart
#    > summary(as.numeric(covid_cantab$datediff.y)/365.25)
#       Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     0.4271  1.2266  1.6619  1.6377  2.0205  3.1814 
#
# 20210310 WF - init

suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(lubridate)
})


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
