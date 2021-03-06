% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.parallel.run.R
\name{model.parallel.run}
\alias{model.parallel.run}
\title{Distribute model runs over multiple cores}
\usage{
model.parallel.run(model, actual.input.matrix, seed_count = 0,
  n_cores = 1)
}
\arguments{
\item{model}{Wrapper function for the model}

\item{actual.input.matrix}{Matrix with parameter combinations to be run}

\item{seed_count}{Origin of random number seed}

\item{n_cores}{Number of cores available for parallel running of the model}
}
\value{
a matrix of model features and the seed of the random number
  generator
}
\description{
Allows faster execution of each wave of simulations
}
