#' MICE with weights in parallel (multiple cores of 1 node)
#'
#' Produces a list with imputations with inverse probability weights.
#'
#' @param data training dataframe
#' @param m number of imputations (1)
#' @param method string (vector)
#' @param predictorMatrix predictorMatrix
#' @param maxit How many cycles through chained equations
#' @param printFlag Do you want stuff printed
#' @param n_cores number of cores available for parallel running of mice
#' @param n_experiments number of experiments (typically a multiple of n.experiments)
#' @return A list with imputations and their attributes
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterApply
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterSetRNGStream

mice.parallel <- function(data,
                          m = 1,
                          method = "norm",
                          predictorMatrix,
                          maxit = 5,
                          printFlag = FALSE,
                          n_cores = 1,
                          n_experiments)
{
  cl <- makeCluster(getOption("cl.cores", n_cores))
  mice.imputation.list <- mice.input.lists <- vector("list", length = n_experiments)

  for (i in 1:n_experiments) {
    #mice.input.list$seed <- seed.init + i
    mice.input.lists[[i]] <- data
  }
  clusterSetRNGStream(cl = cl, 0)

  mice.imputation.list <- clusterApply(cl = cl,
                                       x = mice.input.lists,
                                       fun = mice::mice,
                                       m = m,
                                       method = method,
                                       predictorMatrix = predictorMatrix,
                                       maxit = maxit,
                                       printFlag = printFlag)
  stopCluster(cl)
  return(mice.imputation.list)
}
