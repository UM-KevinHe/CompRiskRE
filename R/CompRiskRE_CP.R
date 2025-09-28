#' Competing Risks with Relative Entropy Integration with a Competing Risk Prior Model
#'
#' This function implements the Relative Entropy (RE) framework for discrete-time competing risk models 
#' with the prior model is a competing risk prior model.
#'
#' @param beta_cor Numeric correlation level between external and local models.
#'   Must be one of \code{1.0, 0.9, 0.7, 0.5, 0.4, 0.3, 0.1, 0}.
#' @param eta Numeric vector of Relative Entropy regularization parameters.
#' @param N_ext Sample size for external data (default: 5000).
#' @param N_loc Sample size for local data (default: 2000).
#' @param N_val Sample size for validation data (default: 200).
#' @param N_test Sample size for test data (default: 5000).
#' @param nCause Number of competing causes (default: 2).
#' @param mu Mean of covariates (default: 1).
#' @param sigma Standard deviation of covariates (default: 0.05).
#' @param seed Random seed for reproducibility.
#'
#' @return A list with two components:
#' \describe{
#'   \item{models}{Fitted models: \code{prior}, \code{local}, \code{joint}, and
#'   \code{KL} (per eta value).}
#'   \item{metrics}{Performance results including validation/test deviance,
#'   C-indices, and AIC across eta values.}
#' }
#'
#' @examples
#' \dontrun{
#' eta <- generate_eta(method = "exponential", n = 30, max_eta = 30)
#' res <- CompRiskRE_CP(beta_cor = 0.9, eta = eta, seed = 2024)
#' res$metrics
#' }
#'
#' @export
CompRiskRE_CP <- function(beta_cor = 0.9,
                          eta,
                          N_ext = 5000,
                          N_loc = 2000,
                          N_val = 200,
                          N_test = 5000,
                          nCause = 2,
                          mu = 1,
                          sigma = 0.05,
                          seed = 123) {
  p <- 4
  Tmax <- 10
  
  # ---- Simulate data ----
  sim_data <- simulate_CPCP(beta_cor, N_ext, N_loc, N_val, N_test,
                            p, Tmax, nCause, mu, sigma, seed)
  

  ext  <- sim_data$external
  loc  <- sim_data$local
  val  <- sim_data$validation
  test <- sim_data$test
  Betat_l <- sim_data$Betat_l; Betav_l <- sim_data$Betav_l
  Betat_p <- sim_data$Betat_p; Betav_p <- sim_data$Betav_p

  
  y_p <- ext$y;  Z_p <- ext$Z
  y_l <- loc$y;  Z_l <- loc$Z
  y_v <- val$y;  Z_v <- val$Z
  y_t <- test$y; Z_t <- test$Z
  
  t_test   <- test$X$t
  ind_test <- test$X$ind
  
  x_test <- data.frame(1, test$X[, 1:p])
  colnames(x_test) <- c("timeInt", paste0("X", 1:p))
  x_test$timeInt <- factor(x_test$timeInt)
  
  # ---- Prior model (competing risk multinomial) ----
  model_prior <- nnet::multinom(y_p ~ . - 1, data = Z_p[, c(2, (1:p) + 3)], maxit = 1000)
  prior_y <- predict(model_prior, Z_l[, c(2, (1:p) + 3)], "prob")
  
  # ---- Local model ----
  model_local <- nnet::multinom(rbind(y_l, y_v) ~ . - 1,
                                data = rbind(Z_l[, c(2, (1:p) + 3)], Z_v[, c(2, (1:p) + 3)]),
                                maxit = 1000)
  
  # ---- Joint model ----
  Z_j <- rbind(Z_p, Z_l, Z_v)
  y_j <- rbind(y_p, y_l, y_v)
  model_joint <- nnet::multinom(y_j ~ . - 1, data = Z_j[, c(2, (1:p) + 3)], maxit = 1000)
  
  
  # ---- True model ---- 
  model_truth1 <- model_local
  start.idx <- which(model_truth1$wts != 0)[1]
  model_truth1$wts[(start.idx):(start.idx + p + Tmax - 1)] <- c(Betat_l[, 1], Betav_l[, 1])
  model_truth1$wts[(start.idx + p + Tmax + 1):length(model_truth1$wts)] <- c(Betat_l[, 2], Betav_l[, 2])
  pred_tr1 <- predict(model_truth1, Z_t[, c(2, (1:p) + 3)], "prob")
  
  KL_vali_rst <- c(); KL_test_rst <- c()
  KL_test_rst_C1 <- c(); KL_test_rst_C2 <- c(); KL_test_rst_ER <- c(); AIC <- c()
  Model_KL <- list()
  
  for (e in eta) {
    y_KL <- (y_l + e * prior_y) / (1 + e)
    model_KL <- nnet::multinom(y_KL ~ . - 1, data = Z_l[, c(2, (1:p) + 3)], maxit = 1000)
    Model_KL[[as.character(e)]] <- model_KL
    
    pred_KL_l <- predict(model_KL, Z_l[, c(2, (1:p) + 3)], "prob")
    AIC <- c(AIC, 2 * nCause * (p + Tmax) / (1 + e) - 2 * sum(y_l * log(pred_KL_l)))
    
    pred_vali_KL <- predict(model_KL, Z_v[, c(2, (1:p) + 3)], "prob")
    KL_vali_rst <- c(KL_vali_rst, -2 * sum(y_v * log(pred_vali_KL)))
    
    pred_test_KL <- predict(model_KL, Z_t[, c(2, (1:p) + 3)], "prob")
    KL_test_rst  <- c(KL_test_rst,  -2 * sum(y_t * log(pred_test_KL)))
    KL_test_rst_ER <- c(KL_test_rst_ER, -2 * sum(pred_tr1 * log(pred_test_KL)))
    
    score_KL_l <- predict(model_KL, x_test, "prob")
    KL_test_rst_C1 <- c(KL_test_rst_C1, CindexCR(t_test, ind_test,
                                                 length(score_KL_l[,2]) - rank(score_KL_l[,2]) + 1, 1))
    KL_test_rst_C2 <- c(KL_test_rst_C2, CindexCR(t_test, ind_test,
                                                 length(score_KL_l[,3]) - rank(score_KL_l[,3]) + 1, 2))
  }
  
  return(list(
    models = list(prior = model_prior, local = model_local,
                  joint = model_joint, KL = Model_KL),
    metrics = list(KL_vali = KL_vali_rst, KL_test = KL_test_rst,
                   KL_test_C1 = KL_test_rst_C1, KL_test_C2 = KL_test_rst_C2,
                   KL_test_ER = KL_test_rst_ER, AIC = AIC)
  ))
}

