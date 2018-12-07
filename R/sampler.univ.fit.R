#' adjusted mice univariate sampler
#'
#' Augmented to also return details of the regression model and residual error sampled.
#'
#' @param data dataframe with missing data in it
#' @param r r
#' @param where points to missing data
#' @param type type
#' @param formula formula
#' @param method method
#' @param yname yname
#' @param k k
#' @param calltype calltype
#' @param user user
#' @param ... ...
#' @return A list with imputations and their attributes, including their inverse probability weights
#'
#' @import mice

sampler.univ.fit <- function(data, r, where, type, formula, method, yname, k, calltype = "type",
          user, ...)
{
  #o#btain.design <- mice:::o#btain.design
  #i#nitialize.chain <- mice:::i#nitialize.chain
  #h#andles.format <- mice:::ha#ndles.format
  #i#s.passive <- mice:::i#s.passive
  #s#ampler.univ <- mice:::s#ampler.univ
  #c#heck.df <- mice:::c#heck.df
  #r#emove.lindep <- mice:::r#emove.lindep


  j <- yname[1L]
  if (calltype == "type") {
    vars <- colnames(data)[type != 0]
    formula <- stats::reformulate(setdiff(vars, j), response = j)
    formula <- stats::update(formula, ". ~ . ")
  }
  if (calltype == "formula") {
    ymove <- setdiff(lhs(formula), j)
    formula <- stats::update(formula, paste(j, " ~ . "))
    if (length(ymove) > 0L)
      formula <- stats::update(formula, paste("~ . + ", paste(ymove,
                                                       collapse = "+")))
  }
  x <- obtain.design(data, formula)
  if (calltype == "type") {
    type <- type[labels(stats::terms(formula))][attr(x, "assign")]
    x <- x[, -1L, drop = FALSE]
    names(type) <- colnames(x)
  }
  if (calltype == "formula") {
    x <- x[, -1L, drop = FALSE]
    type <- rep(1L, length = ncol(x))
    names(type) <- colnames(x)
  }
  y <- data[, j]
  ry <- stats::complete.cases(x, y) & r[, j]
  wy <- stats::complete.cases(x) & where[, j]
  if (all(!wy))
    return(numeric(0))
  cc <- wy[where[, j]]
  if (k == 1L)
    check.df(x, y, ry)
  keep <- remove.lindep(x, y, ry, ...)
  x <- x[, keep, drop = FALSE]
  type <- type[keep]
  if (ncol(x) != length(type))
    stop("Internal error: length(type) != number of predictors")
  f <- paste("mice.impute", method, "fit", sep = ".")
  imputes <- data[wy, j]
  imputes[!cc] <- NA
  args <- c(list(y = y, ry = ry, x = x, wy = wy, type = type),
            user, list(...))

  imputation.list <- do.call(f, args = args)
  imputes[cc] <- imputation.list$imputation
  imputation.list$imputes <- imputes
  imputation.list
}
