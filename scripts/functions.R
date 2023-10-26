#' calculates mean and median with confidence estimates under bootstrapping
#' @param x an sf/data frame/tibble containing at least a grouping variable e.g. species name, and numeric data, such as raw seed yields
#' @param dat_col column containing the numeric attributes for which to compute the summary statistics. 
#' @param name_col column containing the grouping variable, or at least a column to name the output rows
#' @param replicates the number of bootstrap re-sampling events to run
#' @param ... further arguments passed to other functions, notably ?boot, and ?boot.ci, such as the type of confidence interval to calculate (boot.ci, 'type'), whether to run in parallel (boot.ci, 'parallel').
#' @example 
#' spp <- split(iris, f = iris$Species) 
#' lapply(spp, boots, 
#'       dat_col = 'Sepal.Length', name_col = 'Species', replicates = 1000,
#'       type = 'bca',  parallel = "multicore")

boots <- function(x, dat_col, name_col, replicates, ...){
  
  ci_mean_fun <- function(x, i) mean(x[i])
  ci_med_fun <- function(x, i) median(x[i])
  
  s_test_res <- shapiro.test(x[,dat_col]) 
  mean <- boot::boot(data = x[,dat_col], statistic = ci_mean_fun, R = replicates)
  mean_ci <- boot::boot.ci(mean, ...)

  med <- boot::boot(data = x[,dat_col], statistic = ci_med_fun,  R = replicates)
  med_ci <- boot::boot.ci(med, ...)
  
  boot_central <- data.frame(
    taxon = unique(sf::st_drop_geometry(x[,name_col])),
    n = nrow(x), 
    mean = mean[['t0']],
    CI_lwr_mean = mean_ci[['bca']][4],
    CI_upr_mean = mean_ci[['bca']][5],
    median = med[['t0']],
    CI_lwr_med = med_ci[['bca']][4],
    CI_upr_med = med_ci[['bca']][5], 
    Replicates = replicates
  )

  return(boot_central)
  
}


#' fit and select a best variogram for a model
#' 
#' this function is internal to 'modeller', it aims to select a fit variogram, and I am sure does so much better than I do subjectively. 
#' @param x the input data to 'modeller'
#' @param model this is internally computed via initial backwards selection of the maximal models specified to dredge.

varFitter <- function(x, model){
  
  cf = formula('~ x + y')
  
  expone <- nlme:::update.lme(
    model.lm, correlation = nlme::corExp(form = cf, nugget = T))
  gaussian <- nlme:::update.lme(
   model.lm, correlation = nlme::corGaus(form = cf, nugget = T))
  spherical <- nlme:::update.lme(
    model.lm, correlation = nlme::corSpher(form = cf, nugget = T))
  ratio <-nlme:::update.lme(
    model.lm, correlation = nlme::corRatio(form = cf, nugget = T))
  linear <- nlme:::update.lme(
    model.lm, correlation = nlme::corLin(form = cf, nugget = T))
 
  var_term <- model.sel(model.lm, expone, gaussian, spherical, ratio, linear)
  var_form <- row.names(var_term)[1]

  if(grep('expone', var_form)){
    
    final_model <- nlme::lme(
      final_terms, data = x, correlation = nlme::corExp(form = cf, nugget = T), method = 'REML')
    
  } else if(grep('gaussian', var_form)) {
    
    final_model <- nlme::lme(
      final_terms, data = x, correlation = nlme::corGaus(form = cf, nugget = T), method = 'REML')
    
  } else if(grep('spherical', var_form)){
    
    final_model <- nlme::lme(
      final_terms, data = x, correlation = nlme::corSpher(form = cf, nugget = T), method = 'REML')
    
  } else if(grep('ratio', var_form)) {
    
    final_model <- nlme::lme(
      final_terms, data = x, correlation = nlme::corRatio(form = cf, nugget = T), method = 'REML')
    
  } else if(grep('linear', var_form)) { # i think linear is the same as none....
    
    final_model <- nlme::lme(
      final_terms, data = x, correlation = nlme::corLin(form = cf, nugget = T), method = 'REML')
    
  } else { 
    final_model <- nlme::lme(final_terms, data = x, method = 'REML')
    }
  
  return(final_model)
}

predictor <- function(x, y){
  
  # create matrix of possible values. 
  
  ## create quick convex hull to bound Latitude and Longitude with 100 samples. 
  
  bounds <- sf::st_convex_hull(x) |>
    dplyr::mutate(
      Latitude = sf::st_drop_geometry[1],
      Longitude = sf::st_drop_geometry[2]
      )
  
  lats <- 
  longs
  
  
  SPEI_range <- seq(-2.45, 2.5, by = 0.05)
  
  preds <- data.frame(
    loam_prcnt = 1:100, 
    aridity = , 
    bio12 = ,
    bio18 = ,
    ngdd5 = , 
    SPEI6 =  SPEI_range,
    SPEI12 = SPEI_range, 
    SPEI24 = SPEI_range,
    Latitude = , 
    Longitude = 
  )
  
}

length ( seq(-2.4, 2.5, by = 0.1) )


#' fit glm with spatial terms to predict raw seed collection weights
modeller <- function(x, y, outdir, tax_col, ...){
  
  model_terms <- expression(
    SeedMass ~ with(SeedViab) * aridity * MAP * MPWQ * !(SPEI6 && SPEI 12 && SPEI24)
    )
  
  # establish a cluster
  clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
  clust <- try(
    parallel::makeCluster(
      getOption("cl.cores", parallel::detectCores()), type = clusterType)
    )
  clusterExport(clust, x)
  
  models <- MuMIn::pdredge(x, cluster = clust, 
                 subset = model_terms, m.lim = c(0, 4), family = 'poisson')
  
  msAICc <- model.sel(models)
  model <- model.avg(msAICc[ msAICc$delta < 2.0, ], fit = TRUE)
  
  ## hopefully skippable!!!
  
  ## now use these terms to create a spatially explicit model 
  model[['formula']] 
  # or rather if sub is safe, the following retrieves the independent variables. 
  final_terms <- gsub('[(][)]', '', model[['formula']][3])
  
  ## end skippable ??
  
 ## now we assess which, if any, variogram model to include in the final model fit.
  final_model <- varFitter(x, model = model[['formula']]) # pipe directly in... ?

  # save r model object
  
  if(!dir.exists(outdir)){(dir.create(outdir, showWarnings = F))}
  
  model_name <- file.path(outdir, unique(tax_col))
  saveRDS(final_model, model_name)
  
  # predict values from old model onto new. 
  
  
  
  
}


str = model[['formula']]
gsub('', "", str)

gsub('[(][)]', '', str[3])
library(MuMIn)

summary( glm(formula = y ~ X1 + X2 + X3 + X2:X3, data = Cement) )

## ranked with AICc by default
(msAICc <- model.sel(fm1, fm2, fm3))
model <- model.avg( msAICc[ msAICc$delta < 2.0, ], fit = TRUE)

row.names(msAICc)[1]
