#' Distribute model runs over multiple cores
#'
#' Allows faster execution of each wave of simulations
#'
#' @param model Wrapper function for the model
#' @param actual.input.matrix Matrix with parameter combinations to be run
#' @param seed_count Origin of random number seed
#' @param n_cluster Number of cores available for parallel running of the model
#' @return a matrix of model features and the seed of the random number
#'   generator
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel parLapplyLB
#' @importFrom magrittr %>%
#' @export

model.parallel.run <- function(model,
                             actual.input.matrix,
                             seed_count = 0,
                             n_cluster = 4){
  cl <- makeCluster(getOption("cl.cores", n_cluster))
  tab_simul_summarystat = NULL
  list_param <- list(NULL)

  nb_simul <- nrow(actual.input.matrix)

  for (i in 1:nb_simul) {
    param <- c((seed_count + i), actual.input.matrix[i, ])
    list_param[[i]] <- param
  }
  list_simul_summarystat = parLapplyLB(cl, list_param,
                                       model)
  tab_simul_summarystat <- do.call(rbind, list_simul_summarystat)
  stopCluster(cl)
  return(cbind(tab_simul_summarystat, seed_count + 1:nb_simul))
}
