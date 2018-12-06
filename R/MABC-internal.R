### Internal functions, copied from the mice package (mice_3.3.0), as available on https://github.com/cran/mice
# Citation: van Buuren S, Groothuis-Oudshoorn K (2011). “mice: Multivariate Imputation by Chained Equations in R.” Journal of Statistical Software, 45(3), 1-67. https://www.jstatsoft.org/v45/i03/. The mice R package is under License: GPL-2 | GPL-3
# Code was copied as per https://github.com/cran/mice/commit/53f69107bb81f03e98dcdd19e90186043864c670

# The function ma_exists was copied from https://github.com/alexanderrobitzsch/miceadds/blob/master/R/ma_exists.R

check.data <- function(data, method) {
  check.dataform(data)

}

check.dataform <- function(data) {
  if (!(is.matrix(data) || is.data.frame(data)))
    stop("Data should be a matrix or data frame", call. = FALSE)
  if (ncol(data) < 2)
    stop("Data should contain at least two columns", call. = FALSE)
  data <- as.data.frame(data)
  mat <- sapply(data, is.matrix)
  if (any(mat)) stop("Cannot handle columns with class matrix: ",
                     colnames(data)[mat])

  dup <- duplicated(colnames(data))
  if (any(dup)) stop("Duplicate names found: ",
                     paste(colnames(data)[dup], collapse = ", "))

  data
}

check.m <- function(m) {
  m <- m[1L]
  if (!is.numeric(m))
    stop("Argument m not numeric", call. = FALSE)
  m <- floor(m)
  if (m < 1L)
    stop("Number of imputations (m) lower than 1.", call. = FALSE)
  m
}

check.cluster <- function(data, predictorMatrix) {
  # stop if the cluster variable is a factor
  isclassvar <- apply(predictorMatrix == -2, 2, any)
  for (j in colnames(predictorMatrix)) {
    if (isclassvar[j] && lapply(data, is.factor)[[j]])
      stop("Convert cluster variable ", j, " to integer by as.integer()")
  }
  TRUE
}

edit.setup <- function(data, setup,
                       allow.na = FALSE,
                       remove.constant = TRUE,
                       remove.collinear = TRUE,
                       remove_collinear = TRUE,
                       ...) {
  # legacy handling
  if (!remove_collinear) remove.collinear <- FALSE

  # edits the imputation model setup
  # When it detec constant or collinear variables, write in loggedEvents
  # and continues imputation with reduced model

  pred <- setup$predictorMatrix
  meth <- setup$method
  vis <- setup$visitSequence
  post <- setup$post

  # FIXME: this function is not yet adapted to blocks
  if (ncol(pred) != nrow(pred) || length(meth) != nrow(pred)
      || ncol(data) != nrow(pred))
    return(setup)

  varnames <- colnames(data)

  # remove constant variables but leave passive variables untouched
  for (j in seq_len(ncol(data))) {
    if (!is.passive(meth[j])) {
      d.j <- data[, j]
      v <- if (is.character(d.j)) NA else var(as.numeric(d.j), na.rm = TRUE)
      constant <- if (allow.na) {
        if (is.na(v)) FALSE else v < 1000 * .Machine$double.eps
      } else {
        is.na(v) || v < 1000 * .Machine$double.eps
      }
      didlog <- FALSE
      if (constant && any(pred[, j] != 0) && remove.constant) {
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "constant")
        didlog <- TRUE
      }
      if (constant && meth[j] != "" && remove.constant) {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog)
          updateLog(out = out, meth = "constant")
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }

  ## remove collinear variables
  ispredictor <- apply(pred != 0, 2, any)
  if (any(ispredictor)) {
    droplist <- find.collinear(data[, ispredictor, drop = FALSE], ...)
  } else {
    droplist <- NULL
  }
  if (length(droplist) > 0) {
    for (k in seq_along(droplist)) {
      j <- which(varnames %in% droplist[k])
      didlog <- FALSE
      if (any(pred[, j] != 0) && remove.collinear) {
        # remove as predictor
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "collinear")
        didlog <- TRUE
      }
      if (meth[j] != "" && remove.collinear) {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog)
          updateLog(out = out, meth = "collinear")
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }

  if (all(pred == 0L)) stop("nothing left to impute")

  setup$predictorMatrix <- pred
  setup$visitSequence <- vis
  setup$post <- post
  setup$method <- meth
  return(setup)
}

check.method <- function(method, data, where, blocks, defaultMethod) {

  if (is.null(method)) return(make.method(data = data,
                                          where = where,
                                          blocks = blocks,
                                          defaultMethod = defaultMethod))
  nimp <- nimp(where, blocks)

  # expand user's imputation method to all visited columns
  # single string supplied by user (implicit assumption of two columns)
  if (length(method) == 1) {
    if (is.passive(method))
      stop("Cannot have a passive imputation method for every column.")
    method <- rep(method, length(blocks))
    method[nimp == 0] <- ""
  }

  # check the length of the argument
  if (length(method) != length(blocks))
    stop("Length of method differs from number of blocks", call. = FALSE)

  # add names to method
  names(method) <- names(blocks)

  # check whether the requested imputation methods are on the search path
  active.check <- !is.passive(method) & nimp > 0 & method != ""
  passive.check <- is.passive(method) & nimp > 0 & method != ""
  check <- all(active.check) & any(passive.check)
  if (check) {
    fullNames <- rep.int("mice.impute.passive", length(method[passive.check]))
  } else {
    fullNames <- paste("mice.impute", method[active.check], sep = ".")
    if (length(method[active.check]) == 0) fullNames <- character(0)
  }
  notFound <- !vapply(fullNames, exists, logical(1),
                      mode = "function", inherits = TRUE)
  if (any(notFound)) {
    stop(paste("The following functions were not found:",
               paste(fullNames[notFound], collapse = ", ")))
  }

  # type checks on built-in imputation methods
  for (j in names(blocks)) {
    vname <- blocks[[j]]
    y <- data[, vname, drop = FALSE]
    mj <- method[j]
    mlist <- list(m1 = c("logreg", "logreg.boot", "polyreg", "lda", "polr"),
                  m2 = c("norm", "norm.nob", "norm.predict", "norm.boot",
                         "mean", "2l.norm", "2l.pan",
                         "2lonly.pan", "quadratic", "ri"),
                  m3 = c("norm", "norm.nob", "norm.predict", "norm.boot",
                         "mean", "2l.norm", "2l.pan",
                         "2lonly.pan", "quadratic", "logreg", "logreg.boot"))
    cond1 <- sapply(y, is.numeric)
    cond2 <- sapply(y, is.factor) & sapply(y, nlevels) == 2
    cond3 <- sapply(y, is.factor) & sapply(y, nlevels) > 2
    if (any(cond1) && mj %in% mlist$m1)
      warning("Type mismatch for variable(s): ",
              paste(vname[cond1], collapse = ", "),
              "\nImputation method ", mj, " is for categorical data.",
              call. = FALSE)
    if (any(cond2) && mj %in% mlist$m2)
      warning("Type mismatch for variable(s): ",
              paste(vname[cond2], collapse = ", "),
              "\nImputation method ", mj, " is not for factors.",
              call. = FALSE)
    if (any(cond3) && mj %in% mlist$m3)
      warning("Type mismatch for variable(s): ",
              paste(vname[cond3], collapse = ", "),
              "\nImputation method ", mj, " is not for factors with >2 levels.",
              call. = FALSE)
  }
  method[nimp == 0] <- ""
  unlist(method)
}

check.blocks <- function(blocks, data, calltype = "type") {

  data <- check.dataform(data)
  blocks <- name.blocks(blocks)

  # check that all variable names exists in data
  bv <- unique(unlist(blocks))
  notFound <- !bv %in% colnames(data)
  if (any(notFound))
    stop(paste("The following names were not found in `data`:",
               paste(bv[notFound], collapse = ", ")))

  if (length(calltype) == 1L) {
    ct <- rep(calltype, length(blocks))
    names(ct) <- names(blocks)
    attr(blocks, "calltype") <- ct
  }
  else {
    ct <- calltype
    names(ct) <- names(blocks)
    attr(blocks, "calltype") <- ct
  }

  blocks
}

check.predictorMatrix <- function(predictorMatrix,
                                  data,
                                  blocks = NULL) {
  data <- check.dataform(data)

  if (!is.matrix(predictorMatrix))
    stop("predictorMatrix not a matrix", call. = FALSE)
  if (any(dim(predictorMatrix) == 0L))
    stop("predictorMatrix has no rows or columns", call. = FALSE)

  # if we have no blocks, restrict to square predictorMatrix
  if (is.null(blocks)) {
    if (nrow(predictorMatrix) != ncol(predictorMatrix))
      stop(paste("If no blocks are specified, predictorMatrix must",
                 "have same number of rows and columns"),
           call. = FALSE)
    if (is.null(dimnames(predictorMatrix))) {
      if (ncol(predictorMatrix) == ncol(data))
        dimnames(predictorMatrix) <- list(colnames(data), colnames(data))
      else
        stop("Missing row/column names in predictorMatrix", call. = FALSE)
    }
    for (i in row.names(predictorMatrix))
      predictorMatrix[i, grep(i, colnames(predictorMatrix))] <- 0
    return(predictorMatrix)
  }

  # check conforming arguments
  if (nrow(predictorMatrix) > length(blocks))
    stop(paste0("predictorMatrix has more rows (", nrow(predictorMatrix),
                ") than blocks (", length(blocks), ")"),
         call. = FALSE)

  # borrow rownames from blocks if needed
  if (is.null(rownames(predictorMatrix)) &&
      nrow(predictorMatrix) == length(blocks))
    rownames(predictorMatrix) <- names(blocks)
  if (is.null(rownames(predictorMatrix)))
    stop("Unable to set row names of predictorMatrix", call. = FALSE)

  # borrow blocknames from predictorMatrix if needed
  if (is.null(names(blocks)) &&
      nrow(predictorMatrix) == length(blocks))
    names(blocks) <- rownames(predictorMatrix)
  if (is.null(names(blocks)))
    stop("Unable to set names of blocks", call. = FALSE)

  # check existence of row names in blocks
  found <- rownames(predictorMatrix) %in% names(blocks)
  if (!all(found))
    stop("Names not found in blocks: ",
         paste(rownames(predictorMatrix)[!found], collapse = ", "),
         call. = FALSE)

  # borrow colnames from data if needed
  if (is.null(colnames(predictorMatrix)) &&
      ncol(predictorMatrix) == ncol(data))
    colnames(predictorMatrix) <- names(data)
  if (is.null(colnames(predictorMatrix)))
    stop("Unable to set column names of predictorMatrix", call. = FALSE)

  # check existence of variable names on data
  found <- colnames(predictorMatrix) %in% names(data)
  if (!all(found))
    stop("Names not found in data: ",
         paste(colnames(predictorMatrix)[!found], collapse = ", "),
         call. = FALSE)

  list(predictorMatrix = predictorMatrix,
       blocks = blocks)
}

is.passive <- function(string) {
  return("~" == substring(string, 1, 1))
}

check.blots <- function(blots, data, blocks = NULL) {
  data <- check.dataform(data)

  if (is.null(blots)) return(make.blots(data, blocks))

  blots <- as.list(blots)
  for (i in seq_along(blots)) blots[[i]] <- as.list(blots[[i]])

  if (length(blots) == length(blocks) && is.null(names(blots)))
    names(blots) <- names(blocks)
  blots
}

check.df <- function(x, y, ry) {
  # if needed, writes the df warning message to the log
  df <- sum(ry) - ncol(x) - 1
  mess <- paste("df set to 1. # observed cases:", sum(ry), " # predictors:", ncol(x) + 1)
  if (df < 1 && sum(ry) > 0)
    updateLog(out = mess, frame = 4)
}

remove.lindep <- function(x, y, ry, eps = 1e-04, maxcor = 0.99,
                          allow.na = TRUE, frame = 4, ...) {
  # returns a logical vector of length ncol(x)

  if (ncol(x) == 0)
    return(NULL)
  if (eps <= 0)
    stop("\n Argument 'eps' must be positive.")

  # Keep all predictors if we allow imputation of fully missing y
  if (allow.na && sum(ry) == 0) return(rep.int(TRUE, ncol(x)))

  xobs <- x[ry, , drop = FALSE]
  yobs <- as.numeric(y[ry])
  if (var(yobs) < eps) return(rep(FALSE, ncol(xobs)))

  keep <- unlist(apply(xobs, 2, var) > eps)
  keep[is.na(keep)] <- FALSE
  highcor <- suppressWarnings(unlist(apply(xobs, 2, cor, yobs) < maxcor))
  keep <- keep & highcor
  if (all(!keep))
    updateLog(out = "All predictors are constant or have too high correlation.",
              frame = frame)

  # no need to calculate correlations, so return
  k <- sum(keep)
  if (k <= 1L) return(keep)  # at most one TRUE

  # correlation between x's
  cx <- cor(xobs[, keep, drop = FALSE], use = "all.obs")
  eig <- eigen(cx, symmetric = TRUE)
  ncx <- cx
  while (eig$values[k]/eig$values[1] < eps) {
    j <- seq_len(k)[order(abs(eig$vectors[, k]), decreasing = TRUE)[1]]
    keep[keep][j] <- FALSE
    ncx <- cx[keep[keep], keep[keep], drop = FALSE]
    k <- k - 1
    eig <- eigen(ncx)
  }
  if (!all(keep)) {
    out <- paste(dimnames(x)[[2]][!keep], collapse = ", ")
    updateLog(out = out, frame = frame)
  }
  return(keep)
}

check.formulas <- function(formulas, data) {
  formulas <- name.formulas(formulas)
  formulas <- handle.oldstyle.formulas(formulas, data)
  formulas <- lapply(formulas, expand.dots, data)
  # escape if formula is list of two formula's
  if (any(sapply(formulas, is.list))) return(formulas)
  formulas <- lapply(formulas, as.formula)
  formulas
}

updateLog <- function(out = NULL, meth = NULL, frame = 1) {

  # find structures defined a mice() level
  pos_state <- ma_exists("state", frame)$pos
  pos_loggedEvents <- ma_exists("loggedEvents", frame)$pos

  s <- get("state", pos_state)
  r <- get("loggedEvents", pos_loggedEvents)

  rec <- data.frame(it = s$it,
                    im = s$im,
                    dep = s$dep,
                    meth = if(is.null(meth)) s$meth else meth,
                    out = if (is.null(out)) "" else out)

  if (s$log)
    rec <- rbind(r, rec)
  s$log <- TRUE
  assign("state", s, pos = pos_state, inherits = TRUE)
  assign("loggedEvents", rec, pos = pos_loggedEvents, inherits = TRUE)
  return()
}

check.post <- function(post, data) {

  if(is.null(post)) return(make.post(data))

  # check
  if (length(post) != ncol(data))
    stop("length(post) does not match ncol(data)", call. = FALSE)

  # change
  if (is.null(names(post))) names(post) <- colnames(data)

  post
}

check.visitSequence <- function(visitSequence = NULL,
                                data, where = NULL, blocks) {

  if (is.null(names(blocks)) || any(is.na(names(blocks))))
    stop("Missing names in `blocks`.")

  if (is.null(visitSequence)) return(make.visitSequence(data, blocks))

  if (is.null(where)) where <- is.na(data)
  nimp <- nimp(where, blocks)
  if (length(nimp) == 0) visitSequence <- nimp

  if (length(visitSequence) == 1 && is.character(visitSequence)) {
    code <- match.arg(visitSequence, c("roman", "arabic", "monotone",
                                       "revmonotone"))
    visitSequence <- switch(
      code,
      roman = names(blocks)[nimp > 0],
      arabic = rev(names(blocks)[nimp > 0]),
      monotone = names(blocks)[order(nimp)],
      revmonotone = rev(names(blocks)[order(nimp)])
    )
  }

  # legacy handling
  if (is.numeric(visitSequence))
    visitSequence <- colnames(data)[visitSequence]

  # check against names(blocks)
  visitSequence <- visitSequence[is.element(visitSequence, names(blocks))]

  # remove any blocks without missing data
  visitSequence <- names((nimp > 0L)[visitSequence])
  visitSequence
}

check.where <- function(where, data, blocks) {
  if (is.null(where))
    where <- make.where(data, keyword = "missing")

  if (!(is.matrix(where) || is.data.frame(where)))
    if (is.character(where)) return(make.where(data, keyword = where))
  else
    stop("Argument `where` not a matrix or data frame", call. = FALSE)
  if (!all(dim(data) == dim(where)))
    stop("Arguments `data` and `where` not of same size", call. = FALSE)

  where <- as.logical(as.matrix(where))
  if (anyNA(where))
    stop("Argument `where` contains missing values", call. = FALSE)

  where <- matrix(where, nrow = nrow(data), ncol = ncol(data))
  dimnames(where) <- dimnames(data)
  where[, !colnames(where) %in% unlist(blocks)] <- FALSE
  where
}

edit.setup <- function(data, setup,
                       allow.na = FALSE,
                       remove.constant = TRUE,
                       remove.collinear = TRUE,
                       remove_collinear = TRUE,
                       ...) {
  # legacy handling
  if (!remove_collinear) remove.collinear <- FALSE

  # edits the imputation model setup
  # When it detec constant or collinear variables, write in loggedEvents
  # and continues imputation with reduced model

  pred <- setup$predictorMatrix
  meth <- setup$method
  vis <- setup$visitSequence
  post <- setup$post

  # FIXME: this function is not yet adapted to blocks
  if (ncol(pred) != nrow(pred) || length(meth) != nrow(pred)
      || ncol(data) != nrow(pred))
    return(setup)

  varnames <- colnames(data)

  # remove constant variables but leave passive variables untouched
  for (j in seq_len(ncol(data))) {
    if (!is.passive(meth[j])) {
      d.j <- data[, j]
      v <- if (is.character(d.j)) NA else var(as.numeric(d.j), na.rm = TRUE)
      constant <- if (allow.na) {
        if (is.na(v)) FALSE else v < 1000 * .Machine$double.eps
      } else {
        is.na(v) || v < 1000 * .Machine$double.eps
      }
      didlog <- FALSE
      if (constant && any(pred[, j] != 0) && remove.constant) {
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "constant")
        didlog <- TRUE
      }
      if (constant && meth[j] != "" && remove.constant) {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog)
          updateLog(out = out, meth = "constant")
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }

  ## remove collinear variables
  ispredictor <- apply(pred != 0, 2, any)
  if (any(ispredictor)) {
    droplist <- find.collinear(data[, ispredictor, drop = FALSE], ...)
  } else {
    droplist <- NULL
  }
  if (length(droplist) > 0) {
    for (k in seq_along(droplist)) {
      j <- which(varnames %in% droplist[k])
      didlog <- FALSE
      if (any(pred[, j] != 0) && remove.collinear) {
        # remove as predictor
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "collinear")
        didlog <- TRUE
      }
      if (meth[j] != "" && remove.collinear) {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog)
          updateLog(out = out, meth = "collinear")
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }

  if (all(pred == 0L)) stop("nothing left to impute")

  setup$predictorMatrix <- pred
  setup$visitSequence <- vis
  setup$post <- post
  setup$method <- meth
  return(setup)
}

handles.format <- function(fn) {
  # determine whether function fn handles the `format` argument
  f <- get(fn)
  handles.arg(f, "format")
}

handles.arg <- function(f, a = "data") {
  # determine whether function f handles argument a
  if (!is.function(f)) return(FALSE)
  a %in% names(formals(f))
}

initialize.chain <- function(blocks, maxit, m) {
  vars <- unique(unlist(blocks))
  chain <- array(NA, dim = c(length(vars), maxit, m))
  dimnames(chain) <- list(vars,
                          seq_len(maxit),
                          paste("Chain", seq_len(m)))
  chain
}

initialize.imp <- function(data, m, where, blocks, visitSequence,
                           method, nmis, data.init) {
  imp <- vector("list", ncol(data))
  names(imp) <- names(data)
  r <- !is.na(data)
  for (h in visitSequence) {
    for (j in blocks[[h]]) {
      y <- data[, j]
      ry <- r[, j]
      wy <- where[, j]
      imp[[j]] <- as.data.frame(matrix(NA, nrow = sum(wy), ncol = m))
      dimnames(imp[[j]]) <- list(row.names(data)[wy], 1:m)
      if (method[h] != "") {
        for (i in seq_len(m)) {
          if (nmis[j] < nrow(data)) {
            if (is.null(data.init)) {
              imp[[j]][, i] <- mice.impute.sample(y, ry, wy = wy)
            } else {
              imp[[j]][, i] <- data.init[wy, j]
            }
          } else imp[[j]][, i] <- rnorm(nrow(data))
        }
      }
    }
  }
  imp
}

obtain.design <- function(data, formula = ~ .) {

  mf <- model.frame(formula, data = data, na.action = na.pass)
  model.matrix(formula, data = mf)
}

sampler.univ <- function(data, r, where, type, formula, method, yname, k,
                         calltype = "type", user, ...) {
  j <- yname[1L]

  if (calltype == "type") {
    vars <- colnames(data)[type != 0]
    formula <- reformulate(setdiff(vars, j), response = j)
    formula <- update(formula, ". ~ . ")
  }

  if (calltype == "formula") {
    # move terms other than j from lhs to rhs
    ymove <- setdiff(lhs(formula), j)
    formula <- update(formula, paste(j, " ~ . "))
    if (length(ymove) > 0L)
      formula <- update(formula, paste("~ . + ", paste(ymove, collapse = "+")))
  }

  # get the model matrix
  x <- obtain.design(data, formula)

  # expand type vector to model matrix, remove intercept
  if (calltype == "type") {
    type <- type[labels(terms(formula))][attr(x, "assign")]
    x <- x[, -1L, drop = FALSE]
    names(type) <- colnames(x)
  }
  if (calltype == "formula") {
    x <- x[, -1L, drop = FALSE]
    type <- rep(1L, length = ncol(x))
    names(type) <- colnames(x)
  }

  # define y, ry and wy
  y <- data[, j]
  ry <- complete.cases(x, y) & r[, j]
  wy <- complete.cases(x) & where[, j]

  # nothing to impute
  if (all(!wy)) return(numeric(0))

  cc <- wy[where[, j]]
  if (k == 1L) check.df(x, y, ry)

  # remove linear dependencies
  keep <- remove.lindep(x, y, ry, ...)
  x <- x[, keep, drop = FALSE]
  type <- type[keep]
  if (ncol(x) != length(type))
    stop("Internal error: length(type) != number of predictors")

  # here we go
  f <- paste("mice.impute", method, sep = ".")
  imputes <- data[wy, j]
  imputes[!cc] <- NA

  args <- c(list(y = y, ry = ry, x = x, wy = wy, type = type), user, list(...))
  imputes[cc] <- do.call(f, args = args)
  imputes
}

find.collinear <- function(x, threshold = 0.999, ...) {
  nvar <- ncol(x)
  x <- data.matrix(x)
  r <- !is.na(x)
  nr <- apply(r, 2, sum, na.rm = TRUE)
  ord <- order(nr, decreasing = TRUE)
  xo <- x[, ord, drop = FALSE]  ## SvB 24mar2011
  varnames <- dimnames(xo)[[2]]
  z <- suppressWarnings(cor(xo, use = "pairwise.complete.obs"))
  hit <- outer(seq_len(nvar), seq_len(nvar), "<") & (abs(z) >= threshold)
  out <- apply(hit, 2, any, na.rm = TRUE)
  return(varnames[out])
}

# This helper function was copied from
# https://github.com/alexanderrobitzsch/miceadds/blob/master/R/ma_exists.R
ma_exists <- function( x, pos, n_index=1:8)
{
  n_index <- n_index + 1
  is_there <- exists(x, where=pos)
  obj <- NULL
  if (is_there){
    obj <- get(x, pos)
  }
  if (! is_there){
    for (nn in n_index){
      pos <- parent.frame(n=nn)
      is_there <- exists(x, where=pos)
      if (is_there){
        obj <- get(x, pos)
        break
      }
    }
  }
  #--- output
  res <- list( is_there=is_there, obj=obj, pos=pos)
  return(res)
}
