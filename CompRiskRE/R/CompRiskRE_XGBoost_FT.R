#' Competing Risks with Relative Entropy Integration using XGBoost and Composite Failure Time Prior
#'
#' This function implements the Relative Entropy (RE) framework for discrete-time competing risks models
#' where the prior model is specified in a composite failure time formulation and the local model is fit
#' using XGBoost. 
#' @param prior_model A fitted prior multinomial model.
#' @param train_data A \code{data.frame} containing training data.
#' @param test_data A \code{data.frame} containing test data.
#' @param eta Numeric vector of RE regularization parameters.
#' @param xgb_params A list of XGBoost parameters..
#' @param nrounds Integer, number of boosting rounds for XGBoost (default: 100).
#' @param maxiter Integer, maximum number of RE updates in \code{priorFTKL_XGBoost} (default: 1).
#' @param eps Numeric tolerance for convergence in \code{priorFTKL_XGBoost} (default: 1e-6).
#'
#' @return A list with two components:
#' \describe{
#'   \item{models}{A named list of fitted XGBoost models for each supplied \code{eta}.}
#'   \item{PD}{A numeric vector of predictive deviance values on the test data.}
#'   \item{eta}{The vector of \code{eta} values.}
#' }
#'
#' @export
CompRiskRE_XGBoost_FT <- function(prior_model, train_data, test_data, eta, 
                                  xgb_params = list("objective" = "multi:softprob",
                                                    "eval_metric" = "mlogloss",
                                                    "num_class" = 3,
                                                    "max.depth" = 4),
                                  nrounds = 100, maxiter = 1, eps = 1e-6) {
  
  nCause <- 2
  p <- ncol(X_train)
  
  # ---- Prepare train data ----
  t_train <- floor(train_data$time_y / 7) + 1
  ind_train <- train_data$event
  X_train <- as.matrix(train_data[, -c("time_y", "event")])
  X_train[is.na(X_train)] <- 0
  
  x_train_long <- as.data.frame(cbind(X_train, t_train, ind_train))
  Z_train <- dataLong(x_train_long, "t_train", "ind_train")
  
  outcome <- matrix(0, nrow = nrow(Z_train), ncol = nCause + 1)
  for (i in 1:nrow(Z_train)) outcome[i, Z_train$y[i] + 1] <- 1
  
  x_train <- fastDummies::dummy_cols(Z_train[, 2])
  x_train <- cbind(x_train[, -1], Z_train[, c((1:p) + 3)])
  y_train <- rep(c(0, 1, 2), each = nrow(Z_train))
  
  x_aug <- rbind(x_train, x_train, x_train)
  outcome_aug <- rbind(outcome, outcome, outcome)
  Weight <- as.vector(outcome_aug)
  
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x_aug), label = y_train, weight = Weight, nthread = 1)

  prior_pred <- nnet::predict(prior_model, Z_train[, c(2, (1:p) + 3)], type = "prob")
  
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
  
  x_test <- fastDummies::dummy_cols(Z_test[, 2])
  x_test <- cbind(x_test[, -1], Z_test[, c((1:p) + 3)])
  dtest <- xgboost::xgb.DMatrix(data = as.matrix(x_test), label = y_test, nthread = 1)
  
  
  # ---- Fit baseline XGBoost model (eta = 0) ----
  x <- fastDummies::dummy_cols(Z_train[, 2])
  x <- cbind(x[, -1], Z_train[, c((1:p) + 3)])
  y <- rep(c(0, 1, 2), each = nrow(Z_train)) 
  
  x_aug <- rbind(x, x, x)
  outcome_aug <- rbind(outcome, outcome, outcome)
  Weight <- as.vector(outcome_aug)
  
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x_aug), label = y, weight = Weight, nthread = 1)
  
  model0 <- xgboost::xgb.train(params = xgb_params, data = dtrain, nrounds = nrounds, 
                               verbose = 2, watchlist = list(test = dtest, train = dtrain))
  
  results <- list()
  results[["0"]] <- model0
  
  PD <- c(min(model0$evaluation_log$test_mlogloss * test_n2 / test_n * 2))
  for (eta_i in eta) {
    rst <- priorFTKL_XGBoost(model0, outcome, x, y, dtrain, dtest, dtrain, xgb_params, prior_pred, 
                             eta = eta_i, maxiter = maxiter, eps = eps)
    results[[as.character(eta_i)]] <- rst$model_KL
    PD <- c(PD, min(rst$model_KL$evaluation_log$test_mlogloss * test_n2 / test_n * 2))
  }
  
  return(list(models = results, PD = PD, eta = eta))
}



priorFTKL_XGBoost <- function(model, outcome, x, y, data, data_test, data_train, xgb_params, prior_y, eta, maxiter = 100, eps = 1e-4) {
  priorCR_y <- cbind(prior_y[, 1], prior_y[, 2] / 2, prior_y[, 2] / 2)
  dev <- numeric(); dev2 <- numeric()
  
  dev <- c(dev, min(model$evaluation_log$test_mlogloss))
  dev2 <- c(dev2, min(model$evaluation_log$test2_mlogloss))
  
  pred <- xgboost::predict(model, newdata = data_train, iterationrange = c(1, 100))
  dim(pred) <- c(3, nrow(data_train))
  pred_KL_l <- t(pred)
  
  priorCR_y <- cbind(prior_y[, 1],
                     prior_y[, 2] * pred_KL_l[, 2] / (pred_KL_l[, 2] + pred_KL_l[, 3]),
                     prior_y[, 2] * pred_KL_l[, 3] / (pred_KL_l[, 2] + pred_KL_l[, 3]))
  
  for (i in 1:maxiter) {
    out_KL <- (outcome + eta * priorCR_y) / (1 + eta)
    Weight <- as.vector(out_KL)
    
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x), label = y, weight = Weight, nthread = 1)
    model <- xgboost::xgb.train(params = xgb_params, data = dtrain, nrounds = 100, verbose = 0, watchlist = list(test = data_test, test2 = data_train))
    
    dev <- c(dev, min(model$evaluation_log$test_mlogloss))
    dev2 <- c(dev2, min(model$evaluation_log$test2_mlogloss))
    
    pred <- xgboost::predict(model, newdata = data_train, iterationrange = c(1, 100))
    dim(pred) <- c(3, nrow(data_train))
    pred_KL_l <- t(pred)
    
    priorCR_y_old <- priorCR_y
    priorCR_y <- cbind(prior_y[, 1],
                       prior_y[, 2] * pred_KL_l[, 2] / (pred_KL_l[, 2] + pred_KL_l[, 3]),
                       prior_y[, 2] * pred_KL_l[, 3] / (pred_KL_l[, 2] + pred_KL_l[, 3]))
    
    if (sum(abs(priorCR_y - priorCR_y_old)) < eps) return(list(model_KL = model, dev = dev, dev2 = dev2, priorCR_y = priorCR_y))
  }
  return(list(model_KL = model, dev = dev, dev2 = dev2, priorCR_y = priorCR_y))
}
