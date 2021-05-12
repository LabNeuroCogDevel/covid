#!/usr/bin/env Rscript
require(neuroCombat) 
require(lme4)
require(dplyr)
require(ggplot2)
require(tidyr)
require(neuroCombat) 
require(lme4)
require(dplyr)
require(glue)

# treat data as global variable in functions
# this is a large (~180Mb) file
# read and parsed using set_CONNDATA()
CONNDATA <- NULL      # data with e.g. hemi's collapsed and anxitey score merged
CONNDATA_orig <- NULL # raw 180Mb csv file

# DEBUG
debug_me <- function() {
    print("hows the data look")
    set_CONNDATA()
    CONNDATA_info()
    print("EXAMPLE MODELS -- everything put together")
    m <- model_one()
    print(summary(m))

    m2 <- model_one(conn="Putamen_avmPFC")
    print(summary(m2))

    print("HARMONIZE EXAMPLE")
    hmz <- harmonize("Putamen_avmPFC")
    str(hmz)
    print("DEFAULT MODEL")
    print(default_model(hmz))

    print("model harm conn")
    print(summary(model_conn(hmz)))
}

# group by all the the columns we might use later
# but not conn (want to mean that) and not 
# 'srcHemi' and 'lat' (collapsing those vlues)
collapse_hemiside <- function(d){
  d %>% group_by(study, subj, ses_id, age, sex, fd_mean, connName) %>%
    summarise(conn=mean(conn)) %>%
    ungroup
}

# take only the first visit
# from either study
only_1ses <- function(d) {
  d %>%
    group_by(subj) %>%
    mutate(vdate=gsub('.*_', '', ses_id),
           vrank=rank(vdate, ties.method="min")) %>%
    filter(vrank < 2) %>%
    select(-vrank) %>%
    ungroup
}

# add composite anxiety score to data. assume input df has 'subj' id
merge_anxiety <- function(CONNDATA, id="subj"){
  # qualtircs composite score for anxiety
  # merge with all.x to keep all CONNDATA even if there isn't an anxiety composite score
  # maybe later remove those without composite like: ... %>% filter(!is.na(composite))
  #   id covid_conn_long is "subj" and in anxiety_score composite was "ID" but renamed in select above
  anx_composite <- read.csv("data/covid_composite_0321.csv") %>%
      select(subj=ID, composite)

  CONNDATA %>%
    only_1ses %>% 
    merge(anx_composite, all.x=T, by=id)
}
#' CONNDATA takes a while to read. do it once here as a global variable
#' names: 
#  "subj"      "study"     "ses_id"    "age"       "sex"       "fd_mean"  
#  "connName"  "conn"      "vdate"     "composite"
set_CONNDATA <- function() {
  # updates global variable. large ~180Mb file
  if(is.null(CONNDATA_orig)){
      CONNDATA_orig <<- read.csv("data/covidv3_conn_long.csv")
      print(glue("read in {nrow(CONNDATA_orig)} conn data rows, with {length(unique(CONNDATA_orig$subj))} subjs"))
  }

  CONNDATA <<- CONNDATA_orig %>%
      collapse_hemiside %>%
      merge_anxiety %>%
      filter(!is.na(composite))
}

#' compare global CONNDATA and CONNDATA_orig
CONNDATA_info <- function() {
  orig <- list(total_row=nrow(CONNDATA_orig), uniqid=length(unique(CONNDATA_orig$subj)))
  new  <- list(total_row=nrow(CONNDATA), uniqid=length(unique(CONNDATA$subj)))
  print(glue("lost {100*(nrow(CONNDATA_orig)-nrow(CONNDATA))/nrow(CONNDATA_orig)}% rows from orig (long fmt {nrow(CONNDATA)} total)"))
  print(glue("adding composite lost: {orig$uniqid - new$uniqid} lost ids. initially: {orig$uniqid}"))
  print(rbind(data.frame(new)%>%mutate(tbl="w/composite"), data.frame(orig)%>%mutate(tbl="orig")))
}


# narrow large df to just the connection we want to inspect
# @param thisConnName roi-roi pair connectivity name "NAcc_Caudate"
subset_conn <- function(thisConnName){
  CONNDATA %>%
    filter(connName == thisConnName,
           age > 0,
           !is.na(conn))
}

# make 'conn' harnomized value, 'conn.original' has old values
# also subset to zscore(conn) < 3
# from FC via MP
harmonize <- function(connName) {
  this <- subset_conn(connName)
  batch = as.factor(c(this$study))
  conns = unlist(this$conn)
  dat = matrix(conns, 1, length(conns))
  
  gender <- as.factor(c(this$sex))
  sage <- c(this$age)
  mod <- model.matrix(~sage+gender)
  
  data.harmonized <- neuroCombat(dat = dat, batch = batch, mod = mod, eb=FALSE)
  data.harmonized <- unlist(data.harmonized)
  # removes extra values from harmonization to match number of values
  data.harmonized3 <- unname(head(data.harmonized, length(this$age))) 
  
  this$data.harmonized3 <- data.harmonized3
  this$conn.original <- this$conn
  this$conn <- this$data.harmonized3
  
  z <- LNCDR::zscore(this$conn)
  this <- this[which(abs(z) < 3),]
}

# set the default model based on what data we give it
# if we have srcHemi as a column use that and lat
# otherwise assume we've collapsed lat and scrHemi
# without multiple values for each roi-roi pair (connName)
# we'll only have 1 value per subject after harmonizing
# so we will not have a random effecto of subj
default_model <- function(d) {
  if('srcHemi' %in% names(d))
    model <- conn ~ age + sex + fd_mean + srcHemi + lat + (1|study/subj)
  else
    model <- conn ~ composite + age + sex + fd_mean + study
}

model_conn <- function(harmonized, model=NULL) {
  # default to ... the default model
  if(is.null(model)) model <- default_model(harmonized)
  
  # lmer if we have 'subj' in the forumla (indicating random effects)
  # lm if simple linear model
  with_rand_ef <- any(grepl('subj',as.character(model)))
  if(with_rand_ef)
    m <- lmer(data=harmonized, model)
  else
    m <- lm(data=harmonized, model)
}

# see a connections values across age and sex
# in the "raw" conn data
plot_conn <- function(conn) {
  conn_df <- subset_conn(conn)
  ggplot(conn_df) +
    aes(y=conn, x=age, color=sex) +
    geom_point() +
    geom_smooth() +
    cowplot::theme_cowplot() +
    ggtitle(glue("{conn} by age w/sex"))
}

# plot the model
# assuming age and sex are factors we want to look at
# and conn is the dependent variable
plot_model_sexcolor <- function(m) {
  pred <- ggeffects::ggpredict(m, term=c("age", "sex"))
  # ggpredict creates dataframe with outcompe 'predict'
  # and then depending on number of terms
  #  1st: x, 2nd: group, 3rd: facet
  
  # get actual data to look like that
  mdata <- dplyr::rename(m@frame, predicted=conn, x=age, group=sex)
  
  ggplot(pred) +
    aes(y=predicted, x=x, group=group, color=group) +
    geom_point(data=mdata, alpha=.3, aes(shape=study)) + # actuall points
    geom_line(aes(color=group)) +
    geom_ribbon(alpha=.1,
                aes(color=NULL, fill=group,
                    ymin=conf.low, ymax=conf.high)) +
    cowplot::theme_cowplot() +
    ggtitle(glue("{conn} by age w/sex"))
}

# wrap up the two steps and catch errors so we can quickly run on everyone
# even if some models don't work (namely roi paired with it's self)
harmonize_and_model <-function(connName, ...){
  tryCatch(model_conn(harmonize(connName), ...),error=function(e){
    warning(e)
    NULL})
}

# get age pvalue from model
#pval_age <- function(m) summary(m)$coefficients['age','Pr(>|t|)']
pval_composite <- function(m) summary(m)$coefficients['composite','Pr(>|t|)']
#age_issig <- function(m, p=.05) pval_age(m) < p
composite_issig <- function(m, p=.05) pval_composite(m) < p

### START HERE
#' harmonize and lm or lmer for each connection
#' @param model model to use. e.g. conn ~ composite + age + sex + fd_mean + study
#'              if null, get from using default_model() on harmonized output
model_harmonize_all <- function(model=NULL){
  # updating global variable! only do it if it doesn't exist
  if(!exists("CONNDATA") || is.null(CONNDATA)) set_CONNDATA()
  
  # if model is null, default set in model_conn is used
  #  model <- conn ~ age + sex + fd_mean + srcHemi + lat + (1|study/subj)
  all_conns <- unique(CONNDATA$connName)
  all_models <- lapply(all_conns, harmonize_and_model, model=model)
  names(all_models) <- all_conns
  return(all_models)
}

#' debug a model for a given connection
#' @param lm/lmer model default NULL uses default_model() 
#' @param conn   default NULL uses first connName in CONNDATA
#' @global CONNDATA
model_one <- function(model=NULL, conn=NULL) {
    if(is.null(conn)) conn <- first(unique(CONNDATA$connName))
    harmonize_and_model(conn)
}



model_CONNDATA <-function(model=NULL){
  # updating global variable! only do it if it doesn't exist
  if(!exists("CONNDATA") || is.null(CONNDATA)) set_CONNDATA()
  
  # if model is null, default set in model_conn is used
  #  model <- conn ~ age + sex + fd_mean + srcHemi + lat + (1|study/subj)
  all_conns <- unique(CONNDATA$connName)
  all_models <- lapply(all_conns, harmonize_and_model, model=model)
  names(all_models) <- all_conns
  
  # remove models that failed
  # N.B. names(all_models) is no longer the same as all_conns
  all_models <- Filter(function(x) !is.null(x), all_models)
  
  all_pvals <- sapply(all_models, pval_age)
  summary(all_pvals)
  # Min.     1st Qu.    Median      Mean   3rd Qu.      Max.
  #0.0000156 0.2011451 0.4555715 0.4596300 0.7302006 0.9766166
  
  names(all_pvals[all_pvals <= .05])
  # "NAcc_amOFC1"    "NAcc_avmPFC"    "NAcc_sgACC"     "NAcc_rACC"     
  # "Caudate_avmPFC" "Caudate_A32p"   "Caudate_A24cd"  "Caudate_dIa"   
  # "Putamen_avmPFC" "Putamen_A24cd"  "aHPC_avmPFC"    "pHPC_pmOFC1"   
  # "amOFC1_amOFC2"  "avmPFC_rACC"   
  
  
  # could use base::Filter (not dpylyr:filter) to retrict. but we already have
  # all the pavlues, so this is redudant 
  #sig_models <- Filter(function(m) age_issig(m, p=.01), all_models)
  
  sig_models_05 <- all_models[all_pvals <= .05]
  cat(glue("{length(all_conns)} total connections, {length(all_models)} modeled, and {length(sig_models_05)} with sig composite effect (p<.05, uncorrected)\n"))
  
  
#  plot_conn(names(sig_models_05)[1])
#  plot_model_sexcolor(sig_models_05[[1]])
  
  return(all_models)
  return(sig_models_05)
}
