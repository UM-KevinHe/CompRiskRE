#' Competing Risks with Relative Entropy Integration (Composite Prior, XGBoost)
#'
#' This function implements the Relative Entropy (RE) framework for discrete-time competing risks models
#' where the prior model is specified in a composite failure time formulation and both prior/local models
#' are fit using XGBoost.
#'
#' @param prior_model A fitted prior XGBoost model.
#' @param train_data A \code{data.frame} containing training data with columns:
#'   covariates, discrete follow-up time (\code{time_y}), and event indicator (\code{event}).
#' @param test_data A \code{data.frame} containing test data with the same structure as \code{train_data}.
#' @param eta Numeric vector of RE regularization parameters.
#' @param xgb_params A list of XGBoost parameters. Defaults to
#'   \code{list(objective = "multi:softprob", eval_metric = "mlogloss", num_class = 3, max.depth = 4)}.
#' @param nrounds Integer, number of boosting rounds for XGBoost (default: 100).
#'
#' @return A list with components:
#' \describe{
#'   \item{models}{A named list of fitted XGBoost models for \eqn{\eta = 0} and each supplied \code{eta}.}
#'   \item{PD}{A numeric vector of predictive deviance values on the test data.}
#'   \item{prior_PD}{Predictive deviance of the prior model on the test data.}
#'   \item{eta}{The vector of \code{eta} values, including \eqn{0} for the baseline model.}
#' }
#'
#' @export
CompRiskRE_XGBoost_CP <- function(prior_model, train_data, test_data, eta,
                                  xgb_params = list("objective" = "multi:softprob",
                                                    "eval_metric" = "mlogloss",
                                                    "num_class" = 3,
                                                    "max.depth" = 4),
                                  nrounds = 100) {
  nCause <- 2
  p <- ncol(X_train)
  
  # ---- Prepare train data ----
  t_train <- floor(train_data$time_y / 7) + 1
  ind_train <- train_data$event
  X_train <- as.matrix(train_data[, -c("time_y", "event")])
  X_train[is.na(X_train)] <- 0
  
  x_train_long <- as.data.frame(cbind(X_train, t_train, ind_train))
  Z_train <- dataLong(x_train_long, "t_train", "ind_train")
  y_train <- Z_train$y
  
  outcome <- matrix(0, nrow = nrow(Z_train), ncol = nCause + 1)
  for (i in 1:nrow(Z_train)) outcome[i, Z_train$y[i] + 1] <- 1

  x_train <- fastDummies::dummy_cols(Z_train[, 2])
  x_train <- cbind(x_train[, -1], Z_train[, c((1:p) + 3)])
  y <- rep(c(0, 1, 2), each = nrow(Z_train))
  
  x_aug <- rbind(x_train, x_train, x_train)
  outcome_aug <- rbind(outcome, outcome, outcome)
  Weight <- as.vector(outcome_aug)
  
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x_aug), label = y, weight = Weight)
  
  dtrain_tmp <- xgboost::xgb.DMatrix(data = as.matrix(x_train))
  pred_prior <- xgboost::predict(prior_model, newdata = dtrain_tmp)
  dim(pred_prior) <- c(3, nrow(dtrain_tmp))
  pred_prior <- t(pred_prior)
  
  # ---- Prepare test data ----
  t_test <- floor(test_data$time_y / 7) + 1
  ind_test <- test_data$event
  X_test <- as.matrix(test_data[, -c("time_y", "event")])
  X_test[is.na(X_test)] <- 0
  
  x_test_long <- as.data.frame(cbind(X_test, t_test, ind_test))
  Z_test <- dataLong(x_test_long, "t_test", "ind_test")
  y_test <- Z_test$y
  
  outcome_test <- matrix(0, nrow = nrow(Z_test), ncol = nCause + 1)
  for (i in 1:nrow(Z_test)) outcome_test[i, Z_test$y[i] + 1] <- 1
  
  # dummy encoding for test
  x_test <- fastDummies::dummy_cols(Z_test[, 2])
  x_test <- cbind(x_test[, -1], Z_test[, c((1:p) + 3)])
  dtest <- xgboost::xgb.DMatrix(data = as.matrix(x_test), label = y_test)
  
  test_n <- nrow(test_data)
  
  # prior predictive deviance
  pred_prior_test <- xgboost::predict(prior_model, newdata = dtest)
  dim(pred_prior_test) <- c(3, nrow(dtest))
  pred_prior_test <- t(pred_prior_test)
  prior_PD <- -sum(outcome_test * log(pred_prior_test)) / test_n * 2
  
  # ---- Fit baseline XGBoost (eta = 0) ----
  model0 <- xgboost::xgb.train(params = xgb_params, data = dtrain, nrounds = nrounds,
                               verbose = 2, watchlist = list(test = dtest, train = dtrain))
  
  results <- list()
  results[["0"]] <- model0

  # ---- Fit KL models for eta > 0 ----
  for (eta_i in eta) {
    out_KL <- (outcome + eta_i * pred_prior) / (1 + eta_i)
    Weight <- as.vector(out_KL)
    dtrain_eta <- xgboost::xgb.DMatrix(data = as.matrix(x_train), label = y, weight = Weight)
    
    model_eta <- xgboost::xgb.train(params = xgb_params, data = dtrain_eta,
                                    nrounds = nrounds, verbose = 2,
                                    watchlist = list(test = dtest, train = dtrain_eta))
    results[[as.character(eta_i)]] <- model_eta
  }
  
  return(list(models = results, PD = PD, prior_PD = prior_PD, eta = c(0, eta)))
}