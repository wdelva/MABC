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
#' @param start.experiments Matrix of experiments (possibly output from a
#'   previous calibration). Set to 0, start experiments will be drawn uniformly
#'   from the prior distributions.
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
                 start.experiments = 0,
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
  calibration.list <- list() # initiating the list where all the output of MiceABC will be stored
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
      if (identical(start.experiments, 0)) {
        # 2. Initial, naive results, based on Sobol sequences
        range.width <- uls - lls
        ll.mat <- matrix(rep(lls, n.experiments), nrow = n.experiments, byrow = TRUE)
        range.width.mat <- matrix(rep(range.width, n.experiments), nrow = n.experiments, byrow = TRUE)
        sobol.seq.0.1 <- sobol(n = n.experiments, dim = length(lls), init = TRUE, scrambling = 1, seed = 1, normal = FALSE)
        experiments <- ll.mat + sobol.seq.0.1 * range.width.mat
      } else {
        experiments <- start.experiments
      }

    }

    if (multinode == TRUE){
      sim.results.simple <- model.snow.run(model = model,
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

    if (wave==1){ # TRUE for the first wave only
      sim.results.with.design.df <- rbind(sim.results.with.design.df,
                                          new.sim.results.with.design.df)
    } else {
      sim.results.with.design.df <- rbind(dplyr::select(sim.results.with.design.df,
                                                        -contains("RMSD")),
                                          new.sim.results.with.design.df)
    }

    ######## CONTEXT
    # sim.results.with.design.df contains all simulations from previous waves (sim.results.with.design.df)
    # and the simulations from the most recent wave (new.sim.results.with.design.df)
    ########

    # sim.results.with.design.df.median.features <- l1median(dplyr::select(sim.results.with.design.df, contains("y.")))
    ######## CONTEXT
    # sim.results.with.design.df.median.features are the median features of ALL experiments for all waves done so far
    ########

    # 3. Find intermediate features and RMSD.tol for which n.close.to.targets >= min.givetomice
    # targets.diff <- targets.empirical - sim.results.with.design.df.median.features # experim.median.features # First we determine how far the empirical targets are away from the median features of the executed experiments
    #candidate.RMSD.tol <- Inf # Initially, we assume that the RMSD cut-off needs to be infinitely large to have sufficient observations to give to mice.

    ### OLD, SLOW APPROACH
    # # Initiate n.close.to.targets
    # n.close.to.targets <- 0 # This will be overwritten.
    # # candidate.intermediate.features <- targets.empirical # We start with the empirical target features
    # RMSD.tol <- 0 # This will be increased if n.close.to.targets < min.givetomice for this tolerance level
    #
    # while (n.close.to.targets < min.givetomice & RMSD.tol <= RMSD.tol.max){
    #
    #   diff.matrix <- sweep(x = sim.results.with.design.df[ , ((1 + x.offset):(x.offset + length(targets.empirical)))], MARGIN = 2, targets.empirical)
    #   rel.diff.matrix <- sweep(diff.matrix, MARGIN = 2, targets.empirical, FUN = "/")
    #   squared.rel.diff.matrix <- rel.diff.matrix^2
    #   sum.squared.rel.diff <- rowSums(squared.rel.diff.matrix)
    #   RMSD <- sqrt(sum.squared.rel.diff / length(targets.empirical))
    #   n.close.to.targets <- sum(RMSD <= RMSD.tol, na.rm = TRUE)
    #   #n.close.to.targets.mat[(1+steps.intermediate.targets), (1+steps.RMSD.tol)] <- n.close.to.targets
    #   #large.enough.training.df <- n.close.to.targets >= min.givetomice
    #   RMSD.tol <- RMSD.tol + 0.0001  # Increasing RMSD.tol
    # }

    ### NEW, FAST APPROACH
    diff.matrix <- sweep(x = sim.results.with.design.df[ , ((1 + x.offset):(x.offset + length(targets.empirical)))],
                         MARGIN = 2,
                         targets.empirical)
    rel.diff.matrix <- sweep(diff.matrix, MARGIN = 2,
                             targets.empirical, FUN = "/")
    squared.rel.diff.matrix <- rel.diff.matrix^2
    sum.squared.rel.diff <- rowSums(squared.rel.diff.matrix)
    RMSD <- sqrt(sum.squared.rel.diff / length(targets.empirical))
    n.close.to.targets <- min.givetomice
    ###



    sim.results.with.design.df$RMSD <- RMSD
    final.intermediate.features <- targets.empirical

    # 5. Select n.close.to.targets shortest distances
    dist.order <- order(RMSD) # Ordering the squared distances from small to big.
    selected.distances <- dist.order[1:n.close.to.targets]
    sim.results.with.design.df.selected <- sim.results.with.design.df[selected.distances, ]

    calibration.list$new.sim.results.with.design.df[[wave]] <- new.sim.results.with.design.df

    # 5.aaaa Keeping track of medians
    calibration.list$sim.results.with.design.df.median.features[[wave]] <- pcaPP::l1median(dplyr::select(sim.results.with.design.df, contains("y.")))
    # The median of the simulations in the lastest wave
    calibration.list$new.sim.results.with.design.df.median.features[[wave]] <- pcaPP::l1median(dplyr::select(new.sim.results.with.design.df, contains("y.")))
    # The median of the simulations to give to mice
    calibration.list$sim.results.with.design.df.selected.median.features[[wave]] <- pcaPP::l1median(dplyr::select(sim.results.with.design.df.selected, contains("y.")))

    # 5.b. Record highest RMSD value for that the selected experiments
    max.RMSD <- max(sim.results.with.design.df.selected$RMSD)
    calibration.list$max.RMSD[[wave]] <- max.RMSD
    # 5.c. Record n.close.target
    calibration.list$n.close.to.targets[[wave]] <- n.close.to.targets

    # 6. Record selected experiments to give to mice for this wave
    calibration.list$selected.experiments[[wave]] <- sim.results.with.design.df.selected

    mice.test <- list()
    if (max.RMSD <= RMSD.tol.max){
      # 7. Put intermediate features in dataframe format
      final.intermediate.features.df <- as.data.frame(matrix(final.intermediate.features, ncol = length(final.intermediate.features)))
      names(final.intermediate.features.df) <- y.names

      # 8. Prepare dataframe to give to mice: selected experiments plus intermediate features

      ## DEBUGGING:
      # We need to replace full_join with smartbind because there are no NAs if there are matching x. values for the added y. values

      df.give.to.mice <- gtools::smartbind(dplyr::select(sim.results.with.design.df.selected,
                                                         -one_of(c("RMSD", "seed", "wave"))), # adding target to training dataset
                                           final.intermediate.features.df[rep(1:nrow(final.intermediate.features.df),
                                                                              each = 1000 * n.experiments), ])


      #print(df.give.to.mice)
      if (!identical(strict.positive.params, 0)){
        df.give.to.mice[, strict.positive.params] <- log(df.give.to.mice[, strict.positive.params])
      }
      # probability.params <- 14
      if (!identical(probability.params, 0)){
        df.give.to.mice[, probability.params] <- log(df.give.to.mice[, probability.params] / (1 - df.give.to.mice[, probability.params])) # logit transformation
      }

      # 9. Override default predictorMatrix with a sparser matrix
      # Let's think a bit more carefully about which variables should be allowed as input for which input parameters.
      # IN THE FUTURE THIS COULD BE AUTOMATED WITH VARIABLE SELECTION ALGORITHMS.
      # predictorMatrix <- (1 - diag(1, ncol(df.give.to.mice))) # This is the default matrix.

      #### NEW: Using LASSO to create predictorMatrix (and ignoring the one that was given as a function argument)

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

      # print(c(nrow(df.give.to.mice) - n.experiments, "nrows to give to mice"))
      # do imputation
      mice.test <- tryCatch(mice.fit(df.give.to.mice,
                                     m = 1,
                                     method = method,
                                     defaultMethod = method,
                                     predictorMatrix = predictorMatrix.give.to.mice,
                                     maxit = maxit,
                                     printFlag = FALSE,
                                     seed = 0),
                            error = function(mice.err) {
                              return(list())
                            })
    }
    # print(c(length(mice.test), "this is length of mice.test"))
    if (length(mice.test) > 0){

      # 11. Turn mice proposals into a new matrix of experiments

      experiments <-matrix(unlist(mice.test$imp),
                           byrow = FALSE,
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

  print(c("time", calibration.list$secondspassed), quote = FALSE)


  return(calibration.list)
}