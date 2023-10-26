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




library(spdep)
?nb2listw
?lm.LMtests
spdep::lm.LMtests()



#' spatial weights generator
#'
#' Generate spatial weights for a generalized linear model.