helper.sir <- function(index, list_param){ # one element of the list_param is c(id, b, g)
  design <- list_param[[index]]
  library(SimInf)
  u0 <- data.frame(S = 99000,
                   I = 1000,
                   R = 0)
  tspan <- c(0, 50, 75)
  feature.vect <- rep(NA, 2)
  set.seed(design[1])
  model <- SIR(u0 = u0,
               tspan = tspan,
               beta = design[2],
               gamma = design[3])
  result <- run(model, threads = 1)
  feature.vect[1] <- result@U[2, 2]
  feature.vect[2] <- result@U[2, 3]
  return(feature.vect)
}

helper.sir.par <- function(list_param.element){ # one element of the list_param is c(id, b, g)
  design <- list_param.element
  library(SimInf)
  u0 <- data.frame(S = 99000,
                   I = 1000,
                   R = 0)
  tspan <- c(0, 50, 75)
  feature.vect <- rep(NA, 2)
  set.seed(design[1])
  model <- SIR(u0 = u0,
               tspan = tspan,
               beta = design[2],
               gamma = design[3])
  result <- run(model, threads = 1)
  feature.vect[1] <- result@U[2, 2]
  feature.vect[2] <- result@U[2, 3]
  return(feature.vect)
}

