## bootci

ci_mean_fun <- function(x,i) mean(x[i])
ci_med_fun <- function(x,i) median(x[i])

#' calculates mean and median with confidence estimates under bootstrapping
#' @param x an sf/data frame/tibble containing at least a grouping variable e.g. species name, and numeric data, such as raw seed yields
#' @param dat_col column containing the numeric attributes for which to compute the summary statistics. 
#' @param name_col column containing the grouping variable, or at least a column to name the output rows
#' @param ... further arguments passed to other functions, notably ?boot, and ?boot.ci, such as the type of confidence interval to calculate (boot.ci, 'type'), whether to run in parallel (boot.ci, 'parallel'), or number of replicates (boot, 'R')

boots <- function(x, dat_col, name_col, ...){
  
  s_test_res <- shapiro.test(x[,dat_col]) 

  mean <- boot::boot(data = x[,dat_col], statistic = ci_mean_fun,  R = 100, parallel = "multicore")
  return(mean)
#  mean_ci <- boot::boot.ci(mean, type = 'bca')

#  med <- boot::boot(data = x[,dat_col], statistic = ci_med_fun,  R = 100)
#  med_ci <- boot::boot.ci(med, type = 'bca', parallel = "multicore")
  
  infls_plants <- data.frame(
 #   taxon = unique(st_drop_geometry(x[,])),
    n = nrow(x), 
    mean = mean[['t0']],
    CI_lwr_mean = mean_ci[['bca']][4],
    CI_upr_mean = mean_ci[['bca']][5],
    median = median[['t0']],
    CI_lwr_med = median_ci[['bca']][4],
    CI_upr_med = median_ci[['bca']][5] 
  )

  return(infls_plants)
}

boots(iris, dat_col = 'Sepal.Length')

