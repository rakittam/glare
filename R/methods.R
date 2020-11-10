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
logLik.anchorglm <- function(object, newdata=NULL, ...){

  if(!is.null(newdata)) {
    data <- newdata
  } else{
    data <- object$data
  }

  mf <- stats::model.frame(object$formula, data = data)
  Y <- stats::model.response(mf)
  X <- stats::model.matrix(object$formula, data = data)

  b <- object$coefficients
  m <- object$m

  return(object$logLik(b=b,Y=Y,X=X,m=m))
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
coef.anchorglm <- function(object, ...){
  return(object$coefficients)
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
predict.anchorglm <- function(object, newdata=NULL,
                              type = c("link", "response"), ...){

  type <- match.arg(type)

  if(!is.null(newdata)) {
    data <- newdata
  } else{
    data <- object$data
  }

  linkinv <- object$family$linkinv

  b <- object$coefficients
  X <- stats::model.matrix(object$formula, data = data)

  pred <- switch(type,
                 "link"= X%*%b,
                 "response"= linkinv(X%*%b))

  return(pred)
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
residuals.anchorglm <- function(object, newdata=NULL,
                                type = c("deviance", "pearson"), ...){

  type <- match.arg(type)

  if(!is.null(newdata)) {
    data <- newdata
  } else{
    data <- object$data
  }

  res_function <- switch(type,
                         "deviance" = object$devianceRes,
                         "pearson" = object$pearsonRes
  )

  family <- object$family
  linkinv <- family$linkinv

  mf <- stats::model.frame(object$formula, data = data)
  Y <- stats::model.response(mf)
  X <- stats::model.matrix(object$formula, data = data)

  b <- object$coefficients
  m <- object$m

  res <- res_function(b, Y, X, linkinv, m, family, ...)

  return(res)
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
