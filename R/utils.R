#' Concordance Index for Competing Risk Models
#'
#' Computes the concordance index (C-index) for a competing risks setting.
#'
#' @param time A numeric vector of observed survival times (must be positive).
#' @param status An integer vector indicating event status:
#'   \itemize{
#'     \item 0 = censored
#'     \item 1 = cause 1 event
#'     \item 2 = cause 2 event
#'   }
#' @param predicted A numeric vector of predicted event probabilities
#'   corresponding to \code{Cause_int}.
#' @param Cause_int An integer specifying which cause of failure (1 or 2)
#'   to evaluate the concordance index for.
#'
#' @return A numeric scalar giving the concordance index
#'   (value between 0 and 1).
#'
#' @export
CindexCR <- function(time, status, predicted, Cause_int = 1) {
  if (anyNA(time) || anyNA(status) || anyNA(predicted)) {
    stop("The input vectors cannot contain NA")
  }
  if (!all(status %in% c(0, 1, 2))) {
    stop("status must be 0, 1, or 2")
  }
  if (!(Cause_int %in% status)) {
    stop("Invalid input for Cause_int")
  }
  if (min(time) <= 0) {
    stop("Survival time must be positive")
  }
  
  return(CindexCR_cpp(time, status, predicted, Cause_int))
}

#' Generate a Sequence of Tuning Parameters (eta)
#' 
#' @param method Character string specifying the method to generate `eta`.
#'   Options are `"linear"` for a linearly spaced sequence, or `"exponential"`
#'   for an exponentially spaced sequence scaled to `max_eta`. Default is `"exponential"`.
#' @param n Integer, the number of `eta` values to generate. Default is 10.
#' @param max_eta Numeric, the maximum value of `eta` in the sequence. Default is 5.
#'
#' @return Numeric vector of length `n` containing the generated `eta` values.
#'
#' @examples
#' # Generate 10 exponentially spaced eta values up to 5
#' generate_eta(method = "exponential", n = 10, max_eta = 5)
#'
#' # Generate 5 linearly spaced eta values up to 3
#' generate_eta(method = "linear", n = 5, max_eta = 3)
#'
#' @export
generate_eta <- function(method = "exponential", n = 10, max_eta = 5) {
  if (method == "linear") {
    eta_values <- seq(0, max_eta, length.out = n)
  } else if (method == "exponential") {
    eta_values <- exp(seq(log(1), log(100), length.out = n))
    eta_values <- (eta_values - min(eta_values)) / (max(eta_values) - min(eta_values)) * max_eta
  }
  return(eta_values)
}


#' Expand Survival Data into Long Format for Discrete-Time Models
#'
#' @param dataSet A \code{data.frame} containing the survival data.
#' @param timeColumn A character string giving the name of the column with observed survival times.
#' @param censColumn A character string giving the name of the column with the event indicator (0 = censored, 1 = event).
#' @param timeAsFactor Logical; if \code{TRUE}, the time intervals are returned as a factor, otherwise as numeric values.
#' @param remLastInt Logical; if \code{TRUE}, removes observations in the final interval where hazard is always 1.
#' @param aggTimeFormat Logical; if \code{TRUE}, expands to a fixed number of intervals given by \code{lastTheoInt}, regardless of observed time.
#' @param lastTheoInt Integer; the maximum theoretical number of intervals. Required if \code{aggTimeFormat = TRUE}.
#'
#' @return A \code{data.frame} in long format with columns:
#'   \itemize{
#'     \item \code{obj}: subject index
#'     \item \code{timeInt}: discrete time interval
#'     \item \code{y}: binary event indicator at each interval
#'     \item all original variables from \code{dataSet}
#'   }
#'
#' @export
dataLong <- function(dataSet, timeColumn, censColumn,
                     timeAsFactor = TRUE,
                     remLastInt = FALSE,
                     aggTimeFormat = FALSE,
                     lastTheoInt = NULL) {
  
  if (!is.data.frame(dataSet)) stop("dataSet must be a data.frame")
  if (!is.character(timeColumn) || length(timeColumn) != 1) stop("timeColumn must be a character scalar")
  if (!is.character(censColumn) || length(censColumn) != 1) stop("censColumn must be a character scalar")
  if (!timeColumn %in% names(dataSet)) stop("timeColumn not found in dataSet")
  if (!censColumn %in% names(dataSet)) stop("censColumn not found in dataSet")
  if (!all(dataSet[[timeColumn]] == floor(as.numeric(dataSet[[timeColumn]])))) {
    stop("timeColumn must contain only integer values")
  }
  
  c1 <- which(names(dataSet) == timeColumn)
  c2 <- which(names(dataSet) == censColumn)
  
  if (aggTimeFormat) {
    if (is.null(lastTheoInt)) stop("lastTheoInt must be specified if aggTimeFormat=TRUE")
    obj <- rep(seq_len(nrow(dataSet)), each = lastTheoInt)
  } else {
    obj <- rep(seq_len(nrow(dataSet)), dataSet[[c1]])
  }
  
  dataSetLong <- dataSet[obj, , drop = FALSE]
  
  if (aggTimeFormat) {
    timeInt <- if (timeAsFactor) {
      factor(rep(seq_len(lastTheoInt), nrow(dataSet)))
    } else {
      rep(seq_len(lastTheoInt), nrow(dataSet))
    }
    y <- unlist(Map(function(t, e) {
      c(rep(0, t - 1), e, rep(0, lastTheoInt - t))
    }, dataSet[[c1]], dataSet[[c2]]))
  } else {
    timeInt <- unlist(Map(function(t) seq_len(t), dataSet[[c1]]))
    if (timeAsFactor) timeInt <- factor(timeInt)
    y <- unlist(Map(function(t, e) c(rep(0, t - 1), e), dataSet[[c1]], dataSet[[c2]]))
  }
  
  dataSetLong <- cbind(obj = obj, timeInt = timeInt, y = y, dataSetLong)
  
  if (remLastInt) {
    maxInt <- max(as.numeric(as.character(dataSetLong$timeInt)))
    remInd <- which(dataSetLong$y == 1 & dataSetLong$timeInt == maxInt)
    if (length(remInd)) dataSetLong <- dataSetLong[-remInd, ]
  }
  
  dataSetLong
}


#' Simulate event times from a discrete competing risks model
#'
#' Generates event times for multiple competing risks under a discrete-time
#' hazard specification.
#'
#' @param x Covariate matrix (\eqn{N \times p}).
#' @param betat Time-related parameters.
#' @param betav Covariate-related parameters.
#' @param N Number of individuals.
#' @param p Number of covariates.
#' @param Tmax Maximum number of discrete time points.
#' @param nCause Number of competing risks.
#' @param unif Logical; if \code{TRUE}, the last column is sampled uniformly, otherwise from a logistic model.
#'
#' @return A numeric matrix of dimension \eqn{N \times (nCause+1)} where each entry
#'   is the simulated event time (or censoring time) for a given individual and risk.
#'
#' @export
prob_n <- function(x, betat, betav, N, p, Tmax, nCause, unif = TRUE) {
  t <- matrix(0, nrow = N, ncol = nCause + 1)
  denom <- matrix(1, nrow = N, ncol = Tmax)
  
  for (k in 1:nCause) {
    xb <- x %*% betav[, k]
    xbm <- matrix(rep(xb, Tmax), nrow = N, ncol = Tmax)
    bt <- matrix(rep(betat[, k], each = N), nrow = N, ncol = Tmax)
    denom <- denom + exp(bt + xbm)
  }
  
  for (k in 1:nCause) {
    xb <- x %*% betav[, k]
    xbm <- matrix(rep(xb, Tmax), nrow = N, ncol = Tmax)
    bt <- matrix(rep(betat[, k], each = N), nrow = N, ncol = Tmax)
    Fm <- exp(bt + xbm) / denom
    
    for (i in 1:N) {
      for (j in 1:Tmax) {
        temp <- sample(c(0, 1), 1, replace = TRUE, prob = c(1 - Fm[i, j], Fm[i, j]))
        if (temp == 1) {
          t[i, k] <- j
          break
        } else if (j == Tmax) {
          t[i, k] <- Tmax + 1
        }
      }
    }
  }
  
  if (unif) {
    for (i in 1:N) {
      t[i, nCause + 1] <- sample(1:Tmax, 1)
    }
  } else {
    xb <- x %*% betav[, nCause + 1]
    xbm <- matrix(rep(xb, Tmax), nrow = N, ncol = Tmax)
    bt <- matrix(rep(betat[, nCause + 1], each = N), nrow = N, ncol = Tmax)
    Fm <- 1 / (1 + exp(-bt - xbm))
    
    for (i in 1:N) {
      for (j in 1:Tmax) {
        temp <- sample(c(0, 1), 1, replace = TRUE, prob = c(1 - Fm[i, j], Fm[i, j]))
        if (temp == 1) {
          t[i, nCause + 1] <- j
          break
        } else if (j == Tmax) {
          t[i, nCause + 1] <- Tmax
        }
      }
    }
  }
  colnames(t) <- c(paste0("cause", seq_len(nCause)), "censor")
  return(t)
}




