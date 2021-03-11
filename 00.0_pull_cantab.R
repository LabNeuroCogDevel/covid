#!/usr/bin/env Rscript

# read in all MOT and SOC CANTAB sheets (line per event)
# collapse into single sheet per task
#
# will be summerised and combined by 00.1_cantab_merge.R
# end goal is to merge with covid data (00.2_cantab_covid.R)
# 20210310WF - init


suppressPackageStartupMessages({library(dplyr)})

read_one <- function(f){
   # 20210311 - do not want "DESCRIPTION" part of sheet 
   # skip no value line like "","",""... ; end if empty line (line before "DESC...")
   cmd <- glue::glue("perl -lne 'next if /^[\",]+.$/; exit if m/^.$/; print' '{f}'")
   #d <- data.table::fread(cmd=cmd) %>% 
   d <- read.csv(pipe(cmd)) %>% mutate(ld8 = LNCDR::ld8from(f))

   # some MOT CSVs have weird values. column read in as character (?)
   for(probcol in c('hit.latency', 'presentation.number', 'problem.number')){
      if(probcol %in% names(d))
         d[[probcol]] <- as.numeric(d[[probcol]])
   }
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
