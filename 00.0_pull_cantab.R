#!/usr/bin/env Rscript
library(dplyr)

read_one <- function(f){
   d <- read.csv(f) %>% 
      mutate(ld8 = LNCDR::ld8from(f))
   # some MOT CSVs have weird values. column read in as character (?)
   if('hit.latency' %in% names(d)) d$hit.latency <- as.numeric(d$hit.latency)
   return(d)
}

read_all_cantab <- function(glob, limit=Inf) {
    Sys.glob(glob) %>% 
        head(n=limit) %>%
        lapply(read_one) %>%
        bind_rows %>%
        relocate(ld8) %>%
        filter(!is.na(ld8))

}

mot_7t <- read_all_cantab(
    '/Volumes/L/bea_res/Data/Temporary Raw Data/7T/*/*_MOT.csv')
write.csv(mot_7t, 'data/cantab_mot.csv', row.names=F, quote=F)

soc_7t <- read_all_cantab(
    '/Volumes/L/bea_res/Data/Temporary Raw Data/7T/*/*_SOC.csv')
write.csv(soc_7t, 'data/cantab_soc.csv', row.names=F, quote=F)
