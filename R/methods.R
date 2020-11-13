#' log likelihood of a anchor glm object
#'
#' Returns log likelihood of a anchor glm object using train or new data
#'
#' @param object anchor glm object
#' @param newdata used for test data, default is NULL
#' @param ... further arguments passed to or from other methods.
#'
#' @return numeric
#' @export
#' @importFrom stats model.response model.matrix model.frame
logLik.anchorglm <- function(object, newdata=NULL, ...) {

  if (!is.null(newdata)) {
    data <- newdata
  } else {
    data <- object$data
  }

  mf <- model.frame(object$formula, data = data)
  Y <- model.response(mf)
  X <- model.matrix(object$formula, data = data)
  b <- object$coefficients
  m <- object$m

  # Handle different form of input for binomial data
  if (dim(as.matrix(Y))[2] == 2) {
    m <- Y[, 1] + Y[, 2]
    Y <- Y[, 1]
  }

  linkinv <- object$family$linkinv

  object$logLik(b = b, Y = Y, X = X, linkinv = linkinv, m = m)
}

#' coefficients of anchor glm object
#'
#' Returns coefficientes of a anchor glm object
#'
#' @param object anchor glm object
#' @param ... further arguments passed to or from other methods.
#'
#' @return numeric
#' @export
coef.anchorglm <- function(object, ...) {
  object$coefficients
}

#' predictions for anchor glm object
#'
#' Returns prediction on scale of linear predictor or response of a anchor glm object using train or new data
#'
#' @param object anchor glm object
#' @param newdata used for test data, default is NULL
#' @param type the type of prediction required. The default is on the scale of the linear predictors; the alternative "response" is on the scale of the response variable. Thus for a default binomial model the default predictions are of log-odds (probabilities on logit scale) and type = "response" gives the predicted probabilities.
#' @param ... further arguments passed to or from other methods.
#'
#' @return numeric
#' @export
#' @importFrom stats model.matrix
predict.anchorglm <- function(object, newdata=NULL,
                              type = c("link", "response"), ...) {

  type <- match.arg(type)

  if (!is.null(newdata)) {
    data <- newdata
  } else {
    data <- object$data
  }

  linkinv <- object$family$linkinv

  b <- object$coefficients
  X <- model.matrix(object$formula, data = data)

  pred <- switch(type,
                 "link"= X %*% b,
                 "response"= linkinv(X %*% b))

  pred
}

#' residuals of a anchor glm object
#'
#' Returns deviance or pearson residuals of a anchor glm object using train or new data
#'
#' @param object anchor glm object
#' @param newdata used for test data, default is NULL
#' @param type the type of residuals which should be returned. The alternatives are: "deviance" (default) and "pearson". Can be abbreviated.
#' @param ... further arguments passed to or from other methods.
#'
#' @return numeric
#' @export
#' @importFrom stats model.response model.matrix model.frame
residuals.anchorglm <- function(object, newdata=NULL,
                                type = c("deviance", "pearson"), ...) {

  type <- match.arg(type)

  if (!is.null(newdata)) {
    data <- newdata
  } else {
    data <- object$data
  }

  res_function <- switch(type,
                         "deviance" = object$devianceRes,
                         "pearson" = object$pearsonRes
  )

  family <- object$family
  linkinv <- family$linkinv

  mf <- model.frame(object$formula, data = data)
  Y <- model.response(mf)
  X <- model.matrix(object$formula, data = data)
  b <- object$coefficients
  m <- object$m

  # Handle different form of input for binomial data
  if (dim(as.matrix(Y))[2] == 2) {
    m <- Y[, 1] + Y[, 2]
    Y <- Y[, 1]
  }

  res_function(b = b, Y = Y, X = X, linkinv = linkinv, m = m, family = family)
}


# #' Summary of a anchor glm object
# #'
# #' Returns summary of a anchor glm object
# #'
# #' @param object anchor glm object
# #'
# #' @return returns object of class "summary.anchorglm"
# #' @export
# summary.anchorglm <- function(object){
#
#   data <- object$data
#
#
#
#
#
#
#   family <- object$family
#
#   linkinv <- family$linkinv
#   m <- object$m
#   b <- object$optim$par
#   X <- data$X
#   Y <- data$Y
#
#   res <- res_function(b, Y, X, linkinv, m, family, ...)
#
#   return(res)
# }
