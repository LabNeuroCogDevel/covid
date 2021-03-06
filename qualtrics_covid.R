#!/usr/bin/env Rscript

# 20210217WF - move qualtrics to covid directory. data/ instead of txt/
# 20210121WF - combine Kid and Adult covid qualtrics surveys
#  creates  'data/qualtrics_adult-kid.csv' and 'data/covid_battery_name_matches_adult-kid.csv'
#
# INTERACTIVE NOTE:
#  the first time this is run, it is necessary to interactively pair the two suvey questions
#  once this is done, it can be reloaded

# TODO:  currently assumes
#  * only 2 surveys match "Bat.*Covid" and
#  * adult is first, kid is second

suppressPackageStartupMessages({
   library(LNCDR) # 1.4.3 -- 20210121 for qualtRics
   library(jsonlite)
})

# get questionnaires
covids <- LNCDR::qualtrics_surveys('Batt.*Covid|Covid19_Youth')

# regexp created from:
# allsurveys <- qualtRics::all_surveys()
# allsurveys %>% 
#  filter(grepl('covid',name,ignore.case=T)) %>% 
#  select(name, lastModified, isActive) %>%
#  arrange(lastModified)

# surveys may come in any order. (ingially adults was first)
snames <- names(covids)

adults <- covids[[grep('Adult',snames)]]
parent <- covids[[grep('Parent',snames)]]
kids   <- covids[[grep('Covid19_Youth',snames)]]

names(adults) <- LNCDR::qualtrics_labels(adults)
names(parent) <- LNCDR::qualtrics_labels(parent)
names(kids)   <- LNCDR::qualtrics_labels(kids)

# INTERACTIVE: match questions between the two
# but only if we haven't already.
# if match file already exists, reuse it instead of interactive matching
matches_csv <- 'data/covid_battery_name_matches_adult-kid.csv'
if(file.exists(matches_csv)){
    matches <- read.csv(matches_csv)
} else {
    if (!exists('data')) dir.create('data')
    na <- names(adults)
    nk <- names(kids)
    matches_matrix <- LNCDR::interactive_label_match(na, nk, accept_single=T, diffprint=F)
    matches <- as.data.frame(matches_matrix)
    names(matches) <- c('adult','kids')
    write.csv(matches_matrix, matches_csv, row.names=F)
}

# select just the overlapping columns
# and make the names the same
# so we can combine
a <- adults[matches$adult]
k <- kids[matches$kid]
names(k) <- names(a)
all_covid_battery <- rbind(a, k)

# save
write.csv(all_covid_battery, 'data/qualtrics-sharedonly_adult-kid.csv', row.names=F)
# all_covid_battery <- read.csv('data/covid_battery_sharedonly.csv') %>%
#  rename(`External Data Reference`=External.Data.Reference)

# add adult only questions
missing_cols <- ! names(adults) %in% matches[,1]
id_col <- grepl("External.Data.Reference", names(adults))
adult_only <- adults[, missing_cols | id_col ]
all_covid_battery_and_adult <-
    merge(all_covid_battery,
          adult_only,
          by="External Data Reference", all.x=TRUE)
write.csv(all_covid_battery, 'data/qualtrics_adult-kid.csv', row.names=F)
