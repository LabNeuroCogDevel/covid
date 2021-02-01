require(neuroCombat) 

# treat data as global variable in functions
CONNDATA <- read.csv('data/covid_conn_long.csv')


# from FC via MP
subset_conn <- function(connName){
  this <- CONNDATA[CONNDATA$connName == connName,]
  this <- this[this$age >0,]
  this <- this[!is.na(this$conn),]
}

# make 'conn' harnomized value, 'conn.original' has old values
# also subset to zscore(conn) < 3
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

model_conn <- function(harmonized) {
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

model_harmonize <-function(connName){
   tryCatch(model_conn(harmonize(connName)),error=function(e) NULL)
}

# get age pvalue from model
pval_age <- function(m) summary(m)$coefficients['age','Pr(>|t|)']
age_issig <- function(m, p=.05) pval_age(m) < p

### START HERE
model_CONNDATA <-function(){

    all_conns <- unique(CONNDATA$connName)
    all_models <- lapply(all_conns, model_harmonize)
    names(all_models) <- all_conns

    # remove models that failed
    # N.B. names(all_models) is no longer the same as all_conns
    all_models <- Filter(function(x) !is.null(x), all_models)

    plot_conn(names(all_models)[1])
    plot_model_sexcolor(all_models[[1]])

    all_sig_vals <- sapply(all_models, pval_age)
    #sig_models <- Filter(function(m) age_issig(m, p=.01), all_models)
    sig_models_05 <- all_models[all_sig_vals <= .05]
    cat(glue::glue("{length(all_conns)} total connections, {length(all_models)} modeled, and {length(sig_models_01)} (.05) with sig age effect (uncorrected)\n"))

}

