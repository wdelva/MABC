#' MICE with weights in parallel (multiple cores of 1 node)
#'
#' Produces a list with imputations with inverse probability weights.
#'
#' @param data training dataframe
#' @param mice.model wrapper function for mice
#' @param m number of imputations
#' @param n_cores number of cores available for parallel running of mice
#' @param n_experiments number of experiments (typically a multiple of n.experiments)
#' @param method string (vector)
#' @param defaultMethod c("pmm", "logreg", "polyreg", "polr")
#' @param predictorMatrix predictorMatrix
#' @param maxit How many cycles through chained equations
#' @param printFlag Do you want stuff printed
#' @param seed.init seed value to make it reproducible
#' @return A list with imputations and their attributes, including their inverse probability weights
#'
#' @import mice
#' @importFrom parallel makeCluster
#' @importFrom parallel parLapplyLB

mice.fit.parallel <- function(data,
                              mice.model = mice.fit.wrapper,
                              m = 1,
                              n_cores = 1,
                              n_experiments,
                              method = "norm",
                              defaultMethod = "norm",
                              predictorMatrix,
                              maxit = 5,
                              printFlag = TRUE,
                              seed.init = 0)
{
  cl <- makeCluster(getOption("cl.cores", n_cores))
  mice.fit.input.lists <- vector("list", length = n_experiments)
  mice.fit.input.list <- list() # An element of mice.fit.input.lists, used as input for mice.fit

  mice.fit.input.list$data <- data
  mice.fit.input.list$m <- m
  mice.fit.input.list$method <- method
  mice.fit.input.list$predictorMatrix <- predictorMatrix
  mice.fit.input.list$defaultMethod <- defaultMethod
  mice.fit.input.list$maxit <- maxit
  mice.fit.input.list$printFlag <- printFlag

  mice.imputation.list <- NULL

  for (i in 1:n_experiments) {
    mice.fit.input.list$seed <- seed.init + i
    mice.fit.input.lists[[i]] <- mice.fit.input.list
  }
  mice.imputation.list <- parLapplyLB(cl, mice.fit.input.lists,
                                   mice.model)
  stopCluster(cl)
  return(mice.imputation.list)
}
