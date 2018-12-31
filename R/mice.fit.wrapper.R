#' mice.fit wrapper
#'
#' Wraps the mice.fit function for parallel execution.
#'
#' @param mice.fit.input.list list of stuff that mice.fit needs
#' @return A list with imputations and their attributes, including their inverse probability weights
#' @import mice

mice.fit.wrapper <- function(mice.fit.input.list){
  mice.imputation <- mice.fit(data = mice.fit.input.list$data,
                              m = mice.fit.input.list$m,
                              method = mice.fit.input.list$method,
                              defaultMethod = mice.fit.input.list$defaultMethod,
                              predictorMatrix = mice.fit.input.list$predictorMatrix,
                              maxit = mice.fit.input.list$maxit,
                              printFlag = mice.fit.input.list$printFlag,
                              seed = mice.fit.input.list$seed)
  return(mice.imputation)
}
