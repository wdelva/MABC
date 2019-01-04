#' MICE-assisted Approximate Bayesian Calibration
#'
#' Produces a list with multiple waves of proposed input parameter values to
#' match a vector of target features.
#'
#' The \code{model} wrapper function for the simulaton model must have a vector
#' of model input parameter values as its one and only  argument. Furthermore,
#' it must return a vector of model features. These model features are then
#' compared against the target features.
#'
#' @param targets.empirical The vector of target features
#' @param model Wrapper function for the simulation model. See details for a
#'   description of the required format.
#' @param RMSD.tol.max Tolerance for the root mean squared distance between
#'   target features and model output features
#' @param min.givetomice Minimal number of observations in the training dataset
#'   to which MICE is applied
#' @param n.experiments Number of parameter combinations in each wave of model
#' runs
#' @param start.experiments If set to NULL (default), start experiments will be
#'   drawn uniformly from the prior distributions. If a matrix of input
#'   parameter values, possibly output from a previous calibration, models with
#'   these inputs will be run, instead of drawing from the prior distributions.
#'   To resume where a previous calibration ended, you can input
#'   start.experiments as a data.frame with inputs, outputs, seed, wave and RMSD
#'   values.
#' @param lls Vector of lower limits of the prior distribution of input
#'   parameter values
#' @param uls Vector of upper limits of the prior distribution of input
#'   parameter values
#' @param strict.positive.params Vector of indices that indicate which of the
#'   input parameters are strictly positive. Set to zero if there are no such
#'   parameters.
#' @param probability.params Vector of indices that indicate which of the input
#'   parameters are strictly between 0 and 1. Set to zero if there are no such
#'   parameters.
#' @param inside_prior TRUE by default. If FALSE, parameter sampling is not
#'   restricted to the initial ranges of the prior distribution during the
#'   subsequent algorithm steps.
#' @param method Method used by MICE. E.g. "norm" or "rf"
#' @param predictorMatrix Can be "complete", "LASSO", or a user-defined matrix
#'   of indices that indicate which variables are included in the chained
#'   equations in MICE
#' @param maxit The maxit argument used in MICE (number of times that the
#'   chained equations are cycled through)
#' @param maxwaves The maximum number of waves of model runs
#' @param n_cores The number of cores available for parallel model runs. Default = 1, i.e. serial execution of model runs
#' @param multinode TRUE or FALSE (Default). If TRUE, model runs are distributed
#'   over the cores of multiple nodes, using DOsnow and snow as the back-end to
#'   the foreach package. If FALSE and n_cores > 1, model runs are distributed
#'   over the cores of a single node, using the parallel package.
#' @return A list with multiple waves of proposed input parameter values to
#'   match a vector of target features.
#'
#' @import mice
#' @importFrom gsubfn strapplyc
#' @importFrom randtoolbox sobol
#' @importFrom pcaPP l1median
#' @importFrom glmnet cv.glmnet
#' @importFrom mvtnorm dmvnorm
#' @import dplyr
#' @import tidyverse
#' @export

MABC <- function(targets.empirical,
                 model,
                 RMSD.tol.max = 2,
                 min.givetomice = 64,
                 n.experiments = 256,
                 start.experiments = NULL,
                 lls,
                 uls,
                 strict.positive.params,
                 probability.params,
                 inside_prior = TRUE,
                 method = "norm",
                 predictorMatrix = "complete",
                 maxit = 50,
                 maxwaves = 4,
                 n_cores = n_cores,
                 multinode = FALSE){
  # 0. Start the clock
  ptm <- proc.time()
  # initiating the list where all MABC output will be stored
  calibration.list <- vector("list", length = maxwaves)
  wave <- 1 # initiating the loop of waves of simulations (one iteration is one wave)
  max.RMSD <- Inf # initially it is infinitely large, but in later iterations it shrinks
  #sim.results.with.design.df <- NULL # Will be growing with each wave (appending)
  #sim.results.with.design.df.selected <- NULL
  #final.intermediate.features <- NULL

  modelstring <- unlist(paste0(deparse(model), collapse = " "))
  input.vector.length <- max(unique(stats::na.omit(as.numeric(gsubfn::strapplyc(modelstring, "[[](\\d+)[]]", simplify = TRUE))))) - 1 # minus one because the first input parameter is the random seed


  # input.vector.length <- max(unique(na.omit(as.numeric(unlist(strsplit(unlist(paste0(deparse(model), collapse = " ")), "[^0-9]+"))))))
  output.vector.length <- length(targets.empirical)


  x.names <- paste0("x.", seq(1:input.vector.length))
  y.names <- paste0("y.", seq(1:output.vector.length))
  x.offset <- length(x.names)
  sim.results.with.design.df <- data.frame(matrix(vector(), 0, (input.vector.length + output.vector.length),
                                                  dimnames = list(c(), c(x.names, y.names))),
                                           stringsAsFactors = FALSE) # Will be growing with each wave (appending)
  # sim.results.with.design.df.selected I think this does not need to be initiated
  final.intermediate.features <- rep(NA, times = length(targets.empirical))

  # 1. Start loop of waves
  while (wave <= maxwaves & max.RMSD > 0){
    print(c("wave", wave), quote = FALSE)

    if (wave == 1){
      if (identical(start.experiments, NULL)) {
        # 2. Initial, naive results, based on Sobol sequences
        range.width <- uls - lls
        ll.mat <- matrix(rep(lls, n.experiments), nrow = n.experiments, byrow = TRUE)
        range.width.mat <- matrix(rep(range.width, n.experiments), nrow = n.experiments, byrow = TRUE)
        sobol.seq.0.1 <- sobol(n = n.experiments, dim = length(lls), init = TRUE, scrambling = 1, seed = 1, normal = FALSE)
        experiments <- ll.mat + sobol.seq.0.1 * range.width.mat
      } else {
        if (ncol(start.experiments) == input.vector.length) {
          experiments <- start.experiments
        }
      }

    }

    if (exists("experiments", mode = "numeric")){
      if (multinode == TRUE){
        sim.results.simple <- model.mpi.run(model = model,
                                            actual.input.matrix = experiments,
                                            seed_count = 0,
                                            n_cores = n_cores)
      } else {
        sim.results.simple <- model.parallel.run(model = model,
                                                 actual.input.matrix = experiments,
                                                 seed_count = 0,
                                                 n_cores = n_cores)
      }

      new.sim.results.with.design.df <- as.data.frame(cbind(experiments,
                                                            sim.results.simple,
                                                            rep(wave, times = nrow(experiments))))

      names(new.sim.results.with.design.df) <- c(x.names, y.names, "seed", "wave")

      new.sim.results.with.design.complete <- stats::complete.cases(new.sim.results.with.design.df)
      new.sim.results.with.design.df <- dplyr::filter(new.sim.results.with.design.df,
                                                      new.sim.results.with.design.complete)

      if (wave == 1){ # TRUE for the first wave only
        sim.results.with.design.df <- rbind(sim.results.with.design.df,
                                            new.sim.results.with.design.df)
      } else {
        sim.results.with.design.df <- rbind(dplyr::select(sim.results.with.design.df,
                                                          -contains("RMSD")),
                                            new.sim.results.with.design.df)
      }



      diff.matrix <- sweep(x = sim.results.with.design.df[ , ((1 + x.offset):(x.offset + length(targets.empirical)))],
                           MARGIN = 2,
                           targets.empirical)
      rel.diff.matrix <- sweep(diff.matrix, MARGIN = 2,
                               targets.empirical, FUN = "/")
      squared.rel.diff.matrix <- rel.diff.matrix^2
      sum.squared.rel.diff <- rowSums(squared.rel.diff.matrix)
      RMSD <- sqrt(sum.squared.rel.diff / length(targets.empirical))
      sim.results.with.design.df$RMSD <- RMSD

    } else {
      sim.results.with.design.df <- start.experiments
      sim.results.with.design.df$wave <- 1 # everything that happened previously is now part of wave 1
      new.sim.results.with.design.df <- dplyr::select(sim.results.with.design.df,
                                                      -contains("RMSD"))
      RMSD <- sim.results.with.design.df$RMSD
    }

    n.close.to.targets <- min.givetomice
    final.intermediate.features <- targets.empirical

    # 5. Select n.close.to.targets shortest distances
    dist.order <- order(RMSD) # Ordering the squared distances from small to big.
    selected.distances <- dist.order[1:n.close.to.targets]
    sim.results.with.design.df.selected <- sim.results.with.design.df[selected.distances, ]

    calibration.list[[wave]]$new.sim.results.with.design.df <- new.sim.results.with.design.df

    # 5.aaaa Keeping track of medians
    calibration.list[[wave]]$sim.results.with.design.df.median.features <- pcaPP::l1median(dplyr::select(sim.results.with.design.df, contains("y.")))
    # The median of the simulations in the lastest wave
    calibration.list[[wave]]$new.sim.results.with.design.df.median.features <- pcaPP::l1median(dplyr::select(new.sim.results.with.design.df, contains("y.")))
    # The median of the simulations to give to mice
    calibration.list[[wave]]$sim.results.with.design.df.selected.median.features <- pcaPP::l1median(dplyr::select(sim.results.with.design.df.selected, contains("y.")))

    # 5.b. Record highest RMSD value for that the selected experiments
    max.RMSD <- max(sim.results.with.design.df.selected$RMSD)
    calibration.list[[wave]]$max.RMSD <- max.RMSD
    # 5.c. Record n.close.target
    calibration.list[[wave]]$n.close.to.targets <- n.close.to.targets

    # 6. Record selected experiments to give to mice for this wave
    calibration.list[[wave]]$selected.experiments <- sim.results.with.design.df.selected

    mice.test <- list()
    if (max.RMSD <= RMSD.tol.max & wave < maxwaves){
      # 7. Put intermediate features in dataframe format
      final.intermediate.features.df <- as.data.frame(matrix(final.intermediate.features, ncol = length(final.intermediate.features)))
      names(final.intermediate.features.df) <- y.names

      # 8. Prepare dataframe to give to mice: selected experiments plus intermediate features
      df.give.to.mice <- gtools::smartbind(dplyr::select(sim.results.with.design.df.selected,
                                                         -one_of(c("RMSD", "seed", "wave"))), # adding target to training dataset
                                           final.intermediate.features.df)

      if (!identical(strict.positive.params, 0)){
        df.give.to.mice[, strict.positive.params] <- log(df.give.to.mice[, strict.positive.params])
      }
      if (!identical(probability.params, 0)){
        df.give.to.mice[, probability.params] <- log(df.give.to.mice[, probability.params] / (1 - df.give.to.mice[, probability.params])) # logit transformation
      }

      if (is.numeric(predictorMatrix)){
        predictorMatrix.give.to.mice <- predictorMatrix
      }

      if (identical(predictorMatrix, "LASSO")){
        predictorMatrix.LASSO <- diag(0, ncol = ncol(df.give.to.mice), nrow = ncol(df.give.to.mice))
        all.names <- names(df.give.to.mice)

        nrows.training.df <- dplyr::select(sim.results.with.design.df.selected,
                                           -one_of(c("RMSD", "seed", "wave"))) %>% nrow()

        for(y.index in 1:ncol(df.give.to.mice)){
          x4lasso <- as.matrix(df.give.to.mice[1:nrows.training.df, -y.index])
          y4lasso <- as.numeric(df.give.to.mice[1:nrows.training.df, y.index])
          alpha <- 1
          cvfit <- glmnet::cv.glmnet(x = x4lasso,
                                     y = y4lasso,
                                     family = "gaussian",
                                     alpha = alpha,
                                     nlambda = 20)
          remaining.indices <- stats::coef(cvfit, s = "lambda.1se")@i
          nonzero.names <- names(df.give.to.mice[-nrow(df.give.to.mice), -y.index])[remaining.indices] # These are the columns with non-zero coefficients
          col.indices <- all.names %in% nonzero.names
          predictorMatrix.LASSO[y.index, col.indices] <- 1
        }
        predictorMatrix.give.to.mice <- predictorMatrix.LASSO
      }

      if (identical(predictorMatrix, "complete")){
        predictorMatrix.give.to.mice <- (1 - diag(1, ncol(df.give.to.mice)))
      }

      # do imputation
      mice.test <- tryCatch(mice.parallel(df.give.to.mice,
                                          m = 1,
                                          method = method,
                                          predictorMatrix = predictorMatrix.give.to.mice,
                                          maxit = maxit,
                                          printFlag = FALSE,
                                          n_cores = 1,
                                          n_experiments = 2 * n.experiments),
                            error = function(mice.err) {
                              return(list())
                            })
    }
    if (length(mice.test) > 0){

      # 11. Turn mice proposals into a new matrix of experiments
      experiments <- mice.test %>%
        purrr::modify_depth(1, "imp") %>%
        unlist() %>%
        matrix(byrow = TRUE,
               ncol = length(x.names))

      # Before we check the suitability of the new experimental input parameter values, we must backtransform the log values to natural values
      if (!identical(strict.positive.params, 0)){
        experiments[, strict.positive.params] <- exp(experiments[, strict.positive.params])
      }
      # And we must also backtransform the logit-transformed values
      if (!identical(probability.params, 0)){
        experiments[, probability.params] <- exp(experiments[, probability.params]) / (1 + exp(experiments[, probability.params]))
      }

      # Lastly, we must create a new matrix, of the same dimensions as the naive "experiments" matrix,
      # But sampled with replacement using the inverse probability weights.
      experiments.df <- data.frame(experiments)

      # We could add an argument to the function to force the new experiments to respect the boundaries of the prior distributions.
      within.prior.limits <- rep(TRUE, n.experiments)
      if (inside_prior == TRUE){ # experiments.df is a dataframe with n.experiments rows and length(lls) columns
        params.above.lls <- sign(sweep(x = experiments.df, MARGIN = 2, lls)) %>% rowSums()
        params.below.uls <- sign(sweep(x = -experiments.df, MARGIN = 2, -uls)) %>% rowSums()
        within.prior.limits <- params.above.lls %in% length(lls) & params.below.uls %in% length(uls)
        experiments.df <- experiments.df[within.prior.limits, ]
      }
      set.seed(0) # for reproducibility

      # NEW:
      mu.experiments <- colMeans(experiments.df)
      cov.experiments <- stats::cov(experiments.df)
      densities.experiments <- mvtnorm::dmvnorm(experiments.df, mu.experiments, cov.experiments)
      weights.experiments <- 1/densities.experiments
      #

      experiments <- matrix(unlist(dplyr::sample_n(experiments.df,
                                                   size = n.experiments,
                                                   replace = TRUE,
                                                   weight = weights.experiments)),
                            byrow = FALSE,
                            ncol = length(x.names))

      wave <- wave + 1
    } else {
      wave <- maxwaves + 1
    }
  }

  # 15. Target features
  calibration.list$targets.empirical <- targets.empirical

  # 16. Stop clock and return calibration list
  calibration.list$secondspassed <- proc.time() - ptm # Stop the clock

  return(calibration.list)
}
