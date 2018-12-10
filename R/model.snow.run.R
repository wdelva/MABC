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
#' @import foreach
#' @import snow
#' @import doSNOW
#' @export

model.snow.run <- function(model,
                          actual.input.matrix,
                          seed_count = 0,
                          n_cores = 3){

  cl <- snow::makeCluster(n_cores, type="MPI")
  doSNOW::registerDoSNOW(cl)

  nb_simul <- nrow(actual.input.matrix)
  # To avoid the note "no visible binding for global variable ‘irun’" when
  # running R CMD check
  irun <- NULL
  modelfeatures <- foreach(irun = 1:nb_simul,
                           .inorder=TRUE,
                           .combine="rbind") %dopar% {
                             seed <- seed_count + irun
                             model(c(seed, actual.input.matrix[irun, ]))
                           }
  snow::stopCluster(cl)
  modelfeatures.array <- as.array(cbind(modelfeatures,
                                        seed_count + 1:nb_simul))
  dimnames(modelfeatures.array) <- NULL
  return(modelfeatures.array)
}
