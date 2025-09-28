#' Simulate external, local, validation, and test datasets
#'
#' This function generates four datasets (external, local, validation, and test)
#' for discrete-time competing risks models. The data are simulated based on
#' user-specified correlation level \code{beta_cor} between the local and external
#' regression coefficients.
#' 
#' @param beta_cor Numeric correlation level between external and local models.
#'   Must be one of \code{1.0, 0.9, 0.5, 0.1, 0}.
#' @param N_ext Sample size for external data (default: 5000).
#' @param N_loc Sample size for local data (default: 2000).
#' @param N_val Sample size for validation data (default: 200).
#' @param N_test Sample size for test data (default: 5000).
#' @param p Number of covariates (default: 4).
#' @param Tmax Maximum follow-up time (discrete intervals) (default: 10).
#' @param nCause Number of competing causes (default: 2).
#' @param mu Mean of covariates (default: 1).
#' @param sigma Standard deviation of covariates (default: 0.05).
#' @param seed Random seed for reproducibility (default: 123).
#'
#' @return A list containing:
#' \describe{
#'   \item{external}{List with simulated external data: \code{X}, \code{Z}, \code{y}}
#'   \item{local}{List with simulated local data: \code{X}, \code{Z}, \code{y}}
#'   \item{validation}{List with simulated validation data: \code{X}, \code{Z}, \code{y}}
#'   \item{test}{List with simulated test data: \code{X}, \code{Z}, \code{y}}
#'   \item{Betat_l}{Baseline time coefficients for local model}
#'   \item{Betav_l}{Covariate coefficients for local model}
#'   \item{Betat_p}{Baseline time coefficients for external model}
#'   \item{Betav_p}{Covariate coefficients for external model}
#' }
simulate_CPFT <- function(beta_cor = 0.9,
                          N_ext = 5000,
                          N_loc = 2000,
                          N_val = 200,
                          N_test = 5000,
                          p = 4,
                          Tmax = 10,
                          nCause = 2,
                          mu = 1,
                          sigma = 0.05,
                          seed = 123) {
  
  # ---- Local model (fixed) ----
  Betav_l <- cbind(c(-0.05,-0.2,0.05,0.2),
                   c(-0.15,0.1,0.15,-0.1),
                   rep(0, p))
  Betat_l <- cbind(c(-2.3,-2.2,-2.08,-1.95,-1.8,-1.61,-1.39,-1.1,-0.7,0),
                   c(-2.3,-2.2,-2.08,-1.95,-1.8,-1.61,-1.39,-1.1,-0.7,0),
                   rep(0, Tmax))
  
  # ---- External model (Betat_p fixed) ----
  Betat_p <- Betat_l
  if (beta_cor == 1.0) {
    Betav_p <- Betav_l
  } else if (beta_cor == 0.9) {
    Betav_p <- cbind(c(-0.15,-0.1,0.15,0.1),
                     c(-0.25,0.01,0.25,-0.01),
                     rep(0, p))
  } else if (beta_cor == 0.5) {
    Betav_p <- cbind(c(-0.2,-0.05,0.2,0.05),
                     c(-0.3,0.25,0.3,-0.25),
                     rep(0, p))
  } else if (beta_cor == 0.1) {
    Betav_p <- cbind(c(-0.25,-0.01,0.25,0.01),
                     c(-0.35,0.3,0.35,-0.3),
                     rep(0, p))
  } else if (beta_cor == 0) {
    Betav_p <- cbind(c(-0.35,0.1,0.35,-0.1),
                     c(-0.45,0.4,0.45,-0.4),
                     rep(0, p))
  } else {
    stop("beta_cor must be one of 1.0, 0.9, 0.5, 0.1, 0")
  }
  
  # ---- Helper function ----
  generate_data <- function(N, Betat, Betav, sigma_val = sigma) {
    x <- mvtnorm::rmvnorm(N, mean = rep(mu, p), sigma = diag(rep(sigma_val, p)))
    t_mat <- prob_n(x, Betat, Betav, N = N, p = p, Tmax = Tmax, nCause = nCause, unif = TRUE)
    ind <- apply(t_mat, 1, which.min); ind[ind == (nCause + 1)] <- 0
    t <- apply(t_mat, 1, min)
    X <- data.frame(x, t = t, ind = ind)
    Z <- dataLong(X, "t", "ind")
    y <- matrix(0, nrow = nrow(Z), ncol = nCause + 1)
    for (i in 1:nrow(Z)) y[i, Z$y[i] + 1] <- 1
    list(X = X, Z = Z, y = y)
  }
  
  # ---- Generate datasets ----
  set.seed(seed)
  ext  <- generate_data(N_ext,  Betat_p, Betav_p, sigma_val = sigma)
  loc  <- generate_data(N_loc,  Betat_l, Betav_l, sigma_val = sigma)
  val  <- generate_data(N_val,  Betat_l, Betav_l, sigma_val = sigma)
  test <- generate_data(N_test, Betat_l, Betav_l, sigma_val = 1)
  
  return(list(external = ext, local = loc, validation = val, test = test,
              Betat_l = Betat_l, Betav_l = Betav_l,
              Betat_p = Betat_p, Betav_p = Betav_p))
}



simulate_CPCP <- function(beta_cprr = 0.9,
                          N_ext = 5000,
                          N_loc = 2000,
                          N_val = 200,
                          N_test = 5000,
                          p = 4,
                          Tmax = 10,
                          nCause = 2,
                          mu = 1,
                          sigma = 0.05,
                          seed = 123) {
  
  # ---- Local model (fixed) ----
  Betav_l <- cbind(c(-0.25,-1,0.25,1),
                   c(0.25,1,-0.25,-1),
                   rep(0, p))
  Betat_l <- cbind(rep(-6, 10), rep(-6, 10), rep(0, Tmax))
  Betat_l[1:9, 1:2] <- -c(6, 5.5, 5.5, 5.5, 5.5, 5, 5, 5, 5)
  
  # direction matrix for perturbation
  delta <- cbind(c(-0.25, 0.5, 0.25, -0.5),
                 c(0.25, -0.5, -0.25, 0.5),
                 rep(0, p))
  
  # ---- External model ----
  Betat_p <- Betat_l
  if (beta_cprr == 1) {
    Betav_p <- Betav_l
  } else if (beta_cprr == 0.9) {
    Betav_p <- Betav_l + 0.5 * delta
  } else if (beta_cprr == 0.7) {
    Betav_p <- Betav_l + 0.7 * delta
  } else if (beta_cprr == 0.5) {
    Betav_p <- Betav_l + 1.0 * delta
  } else if (beta_cprr == 0.4) {
    Betav_p <- Betav_l + 1.1 * delta
  } else if (beta_cprr == 0.3) {
    Betav_p <- Betav_l + 1.2 * delta
  } else if (beta_cprr == 0.1) {
    Betav_p <- Betav_l + 1.3 * delta
  } else if (beta_cprr == 0) {
    Betav_p <- Betav_l + 1.5 * delta
  } else {
    stop("beta_cprr must be one of 1, 0.9, 0.7, 0.5, 0.4, 0.3, 0.1, 0")
  }
  
  # ---- Helper function ----
  generate_data <- function(N, Betat, Betav, sigma_val = sigma) {
    x <- mvtnorm::rmvnorm(N, mean = rep(mu, p), sigma = diag(rep(sigma_val, p)))
    t_mat <- prob_n(x, Betat, Betav, N = N, p = p, Tmax = Tmax, nCause = nCause, unif = TRUE)
    ind <- apply(t_mat, 1, which.min); ind[ind == (nCause + 1)] <- 0
    t <- apply(t_mat, 1, min)
    X <- data.frame(x, t = t, ind = ind)
    Z <- dataLong(X, "t", "ind")
    y <- matrix(0, nrow = nrow(Z), ncol = nCause + 1)
    for (i in 1:nrow(Z)) y[i, Z$y[i] + 1] <- 1
    list(X = X, Z = Z, y = y)
  }
  
  # ---- Generate datasets ----
  set.seed(seed)
  ext  <- generate_data(N_ext,  Betat_p, Betav_p, sigma_val = sigma)
  loc  <- generate_data(N_loc,  Betat_l, Betav_l, sigma_val = sigma)
  val  <- generate_data(N_val,  Betat_l, Betav_l, sigma_val = sigma)
  test <- generate_data(N_test, Betat_l, Betav_l, sigma_val = 1)
  
  return(list(external = ext, local = loc, validation = val, test = test,
              Betat_l = Betat_l, Betav_l = Betav_l,
              Betat_p = Betat_p, Betav_p = Betav_p))
}