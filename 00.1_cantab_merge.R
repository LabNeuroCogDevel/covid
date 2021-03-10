#!/usr/bin/env Rscript

#
# summerise and combined CANTAB tasks (MOT and SOC) for two studies (PET and 7T)
# depends on 
#   7T:  individual sheets combined by 00.0_pull_cantab.R
#   PET: and on exported summary datasheet (from cantab gui) 
# 20210310WF

library(dplyr)
library(tidyr)
library(lubridate)

yes_to_1 <- function(x) ifelse(x=="yes", 1, 0) # currently no longer used

add_sesnum <- function(d, idcol='ID', dcol='vdate', sescol='session')
    d %>%
        group_by(.data[[idcol]]) %>%
        mutate({{sescol}} := rank(.data[[dcol]]))

# Data
mot_7t <- read.csv('data/cantab_mot.csv')
soc_7t <- read.csv('data/cantab_soc.csv')
pet    <- read.csv('data/cantabsummarydatasheet_pet_20200130.csv')

# merge
mot <- mot_7t %>%
    # mutate(hit=yes_to_1(hit)) %>%
    group_by(ld8) %>%
    summarise(mot_lat=mean(hit.latency, na.rm=T))

soc <- soc_7t %>%
    # mutate(solved.in.min=yes_to_1(solved.in.min)) %>%
    group_by(ld8) %>%
    summarise(soc_num.moves = sum(solved.in.min=="yes", na.rm=T))

motsoc_7t <- merge(mot, soc, by='ld8', all=T) %>%
    separate(ld8, c('ID', 'vdate')) %>%
    mutate(vdate=ymd(vdate)) %>%
    # add session/visit number
    add_sesnum()

motsoc_pet <- pet %>%
    select(ID=Subject.ID,
           vdate=Session.start.time,
           mot_lat=MOT.Mean.latency,
           soc_num.moves=SOC.Problems.solved.in.minimum.moves) %>%
    # remove time part and coearse into Date type
    mutate(vdate=gsub(' \\d+:.*', '', vdate) %>% mdy()) %>%
    add_sesnum()

pet_7t <- rbind(motsoc_pet %>% mutate(study='pet'),
                motsoc_7t %>% mutate(study='7T')) %>%
    # another session column 'anyvisit_num' that counts 7T and PET together
    add_sesnum(sescol='anyvisit_num') %>%
    arrange(ID, vdate)

write.csv(pet_7t, 'data/cantab_pet7T_motsoc.csv', row.names=F, quote=F)
