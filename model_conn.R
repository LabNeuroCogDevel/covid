require(neuroCombat) 
require(lme4)

# treat data as global variable in functions
# this is a large (~180Mb) file
CONNDATA <- read.csv('data/covid_conn_long.csv')


# narrow to large DF to just the connection we want to inspect
# connName like "NAcc_Caudate"
subset_conn <- function(connName){
  this <- CONNDATA[CONNDATA$connName == connName,]
  this <- this[this$age >0,]
  this <- this[!is.na(this$conn),]
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

model_conn <- function(harmonized, model=NULL) {
  # default model
  if(is.null(model))
    model <- conn ~ age + sex + fd_mean + srcHemi + lat + (1|study/subj)

  m <- lmer(data=harmonized, model)
}

plot_conn <- function(conn) {
  conn_df <- subset_conn(conn)
  ggplot(conn_df) +
      aes(y=conn, x=age, color=sex) +
      geom_point() +
      geom_smooth() +
      cowplot::theme_cowplot() +
      ggtitle(glue::glue("{conn} by age w/sex"))
}
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
      ggtitle(glue::glue("{conn} by age w/sex"))
}

harmonize_and_model <-function(connName, ...){
   tryCatch(model_conn(harmonize(connName), ...),error=function(e) NULL)
}

# get age pvalue from model
pval_age <- function(m) summary(m)$coefficients['age','Pr(>|t|)']
age_issig <- function(m, p=.05) pval_age(m) < p

### START HERE
model_CONNDATA <-function(model=NULL){
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
    cat(glue::glue("{length(all_conns)} total connections, {length(all_models)} modeled, and {length(sig_models_05)} with sig age effect (p<.05, uncorrected)\n"))

    
    plot_conn(names(sig_models_05)[1])
    plot_model_sexcolor(sig_models_05[[1]])

    summary(all_models$NAcc_amOFC1)

    return(sig_models_05)
}

