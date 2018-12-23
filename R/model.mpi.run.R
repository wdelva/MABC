#' Distribute model runs over multiple nodes
#'
#' Allows faster execution of each wave of simulations
#'
#' @param model Wrapper function for the model
#' @param actual.input.matrix Matrix with parameter combinations to be run
#' @param seed_count Origin of random number seed
#' @param n_cores Number of slave workers available for parallel running of the model
#' @return a matrix of model features and the seed of the random number
#'   generator
#' @import Rmpi
#' @export

model.mpi.run <- function(model,
                          actual.input.matrix,
                          seed_count = 0,
                          n_cores = 3){
  nb_simul <- nrow(actual.input.matrix)
  list_param <- list(NULL)
  for (i in 1:nb_simul) {
    param <- c((seed_count + i), actual.input.matrix[i, ])
    list_param[[i]] <- param
  }
  mpi.spawn.Rslaves(nslaves = n_cores)
  mpi.bcast.cmd(library(SimInf))
  mpi.bcast.Robj2slave(model)
  modelfeatures.list <- mpi.iapplyLB(1:nb_simul,
                                     model,
                                     list_param = list_param)
  modelfeatures <- do.call(rbind, modelfeatures.list)
  mpi.close.Rslaves()
  mpi.exit()
  return(modelfeatures)
}
