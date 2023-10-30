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
#' @export
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

library(sp)
library(gstat)
data(meuse)
coordinates(meuse) = ~x+y
vgm1 <- variogram(log(zinc)~1, meuse)
fit.variogram(vgm1, vgm(1, "Sph", 300, 1))

#' predict the results from a fit model onto a matrix
#' 
#' This function uses spatial and real time climate values to predict the results of a
#' fit model into a new matrix for immediate assessment of raw seed collection weights. 
#' This function is called within 'modeller'.
#' @param x the data frame from which the model was trained, this used to prevent extreme extrapolation.
#' @param model a model to predict
#' @param vals the number of values per variable for the prediction, defaults to 100 per variable
predictor <- function(x, model, vals){
  
  if(missing(vals)){vals <- 100}
  # create matrix of possible values. 
  
  ## create quick convex hull to bound Latitude and Longitude.  
  bbox <- sf::st_transform(x, 5070) |>
    sf::st_bbox() 

  SPEI_range <- seq(-2.45, 2.5, length.out = vals)
  
  preds <- data.frame(
    viable_prcnt = seq(min(x$viable_prcnt), max = max(x$viable_prcnt), length.outs = vals),
    loam_prcnt = seq(min(x$loam), max(x$loam), length.out = vals), 
    bio1 = seq(min(x$bio1), max(x$bio1), length.out = vals), 
    bio12 = seq(min(x$bio12), max(x$bio12), length.out = vals),
    bio18 = seq(min(x$bio18), max(x$bio18), length.out = vals),
    ngd5 = seq(min(x$ngd5), max(x$ngd5), length.out = vals), 
    SPEI6 =  SPEI_range, # create full stack, model will pull relevant var as
    SPEI12 = SPEI_range,  # required.
    SPEI24 = SPEI_range,
    Latitude = seq(from = bbob['ymin'], to = bbox['ymax'], length.out = vals), 
    Longitude = seq(from = bbob['xmin'], to = bbox['xmax'], length.out = vals)
  ) 
  
}


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
  # subset to only models with a single SPEI value ??? - or replace skippable steps below
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
  if(!dir.exists(file.path(outdir, 'models'))){
    (dir.create(file.path(outdir, 'models'), showWarnings = F))
    }
  
  model_name <- file.path(outdir, unique(tax_col))
  saveRDS(final_model, model_name)
  
  # predict values for a lookup table 
  if(!dir.exists(file.path(outdir, 'predictions'))){
    (dir.create(file.path(outdir, 'predictions'), showWarnings = F))
  }
  preds <- predictor(x, model = final_model, ...)
  
  # save lookup table
  write.csv(preds, file =
              file.path(outdir, 'predictions', paste0(unique(tax_col), '.csv')),
            append = F,  row.names = F)
  message('Predictions for ', unique(tax_col), ' written to disk.')
  
}


#' find the most similar conditions from observed data to a models predictions
#' 
#' The input is a data frame with columns showing the observed values for each variable in the prediction model. This function will load a saved lookup table from disk to determine the most similar prediction. 
#' @return a two row dataframe, one row containing the input information, and the other the matched prediction information
#' @param x a data frame (or data frame/sf/tibble), with all rows identical to the names in the prediction stack
#' @param path a path to a directory containing written out prediction tables, as created by the 'modeller' function.
#' @example 
#' data.frame(viable_prcnt = 1:100, 
#' )
#' @export
most_similar <- function(x, path){
  
  columns <- c(
    'viable_prcnt', 'loam_prcnt', 'bio1', 'bio12', 'bio18', 'ngd5', 
    'SPEI6', 'SPE12', 'SPEI24', 'Latitude', 'Longitude')
  if(missing(path)){path <- '.'} # use current working directory if none supplied.
  
  # check for whether user specified full path to the sub directory or not
  if(grep('predictions$', path)){path} else {file.path(path, 'predictions')}
  
  observed <- sf::st_drop_geometry(ob)
  
  # identify species under consideration
  tax <- observed[,'Taxon']
  
  # load lookup tables
  f <- list.files(path, pattern = '.csv.')
  predicted <- read.csv(file.path( path, f[ grepl(tax, f, ignore.case = TRUE)]))
  
  # identify the most similar row in the data set - rescale variables to avoid undue
  # influence of some with large variation, e.g. lat/long
  predicted <- which.min(
    dist(
      scale(
        rbind(observed, predictions)
      )
    ) 
    [1:nrow(predictions)]
  )
  
  # show variation between the prediction and observed measurements
  p <- dplyr::select(predicted, any_of(columns))
  o <- dplyr::select(observed, any_of(columns))
  change <- ((p / o) / o) * 100 # just put a percent on it. 
  
  output <- dplyr::bind_rows(observed, predicted, change) |> 
    dplyr::mutate(Taxon = tax, .before = 1) |> 
    dplyr::mutate(Metric = c('Observed', 'Predicted', 'Difference')) |> 
    dplyr::relocate(.cols = any_of(columns), .after = Metric) |>
    dplyr::select(-Latitude, -Longitude)
    
  return(output)
  
}

predictions <- data.frame(
  viable_prcnt = seq(from = 1, to = 100, length.out = 500), 
  loam_prcnt = seq(from = 25, to = 45, length.out = 500), 
  bio1 = seq(from = 7, to = 20, length.out = 500),
  bio12 = seq(from = 7, to = 20, length.out = 500),
  bio18 = seq(from = 5, to = 30, length.out = 500),
  ngd5 = seq(from = 100, to = 225, length.out = 500), 
  SPEI6 = seq(from = -2, to = 2, length.out = 500),
  SPEI12 = seq(from = -2, to = 2, length.out = 500),
  SPEI24 = seq(from = -2, to = 2, length.out = 500),
  Latitude = seq(from = 1, to = 1000, length.out = 500), 
  Longitude = seq(from = 1, to = 1000, length.out = 500)
)


observed <- predicted[124,]



