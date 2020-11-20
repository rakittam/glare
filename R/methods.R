#' log likelihood of a anchor glm object
#'
#' Returns log likelihood of a anchor glm object using train or new data
#'
#' @param object anchor glm object
#' @param newdata used for test data, default is NULL
#' @param parameter parameter that should be used, default takes the parameter
#'  from the glare objective.
#' @param ... further arguments passed to or from other methods.
#'
#' @return numeric
#' @export
#' @importFrom stats model.response model.matrix model.frame
logLik.glare <- function(object, newdata = NULL, parameter = NULL, ...) {

  if (!is.null(newdata)) {
    data <- newdata
  } else {
    data <- object$data
  }

  mf <- model.frame(object$formula, data = data)
  Y <- model.response(mf)
  X <- model.matrix(object$formula, data = data)

  if (!is.null(parameter)) {
    b <- parameter
  } else {
    b <- object$coefficients
  }

  # Handle different form of input for binomial data
  if (dim(as.matrix(Y))[2] == 2) {
    m <- Y[, 1] + Y[, 2]
    Y <- Y[, 1]
  } else {
    if (!is.null(newdata)) {
      m <- newdata$m
    } else {
      m <- object$m
    }
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
coef.glare <- function(object, ...) {
  object$coefficients
}

#' predictions for anchor glm object
#'
#' Returns prediction on scale of linear predictor or response of a anchor glm
#'  object using train or new data
#'
#' @param object anchor glm object
#' @param newdata used for test data, default is NULL
#' @param parameter parameter that should be used, default takes the parameter
#'  from the glare objective.
#' @param type the type of prediction required. The default is on the scale of
#'  the linear predictors; the alternative "response" is on the scale of the
#'  response variable. Thus for a default binomial model the default predictions
#'  are of log-odds (probabilities on logit scale) and type = "response" gives
#'  the predicted probabilities.
#' @param ... further arguments passed to or from other methods.
#'
#' @return numeric
#' @export
#' @importFrom stats model.matrix
predict.glare <- function(object, newdata = NULL, parameter = NULL,
                              type = c("link", "response"), ...) {

  type <- match.arg(type)

  if (!is.null(newdata)) {
    data <- newdata
  } else {
    data <- object$data
  }

  linkinv <- object$family$linkinv

  if (!is.null(parameter)) {
    b <- parameter
  } else {
    b <- object$coefficients
  }

  X <- model.matrix(object$formula, data = data)

  pred <- switch(type,
                 "link"= X %*% b,
                 "response"= linkinv(X %*% b))

  pred
}

#' residuals of a anchor glm object
#'
#' Returns deviance or pearson residuals of a anchor glm object using train or
#'  new data
#'
#' @param object anchor glm object
#' @param newdata used for test data, default is NULL
#' @param parameter parameter that should be used, default takes the parameter
#'  from the glare objective.
#' @param type the type of residuals which should be returned. The alternatives
#'  are: "deviance" (default) and "pearson". Can be abbreviated.
#' @param ... further arguments passed to or from other methods.
#'
#' @return numeric
#' @export
#' @importFrom stats model.response model.matrix model.frame
residuals.glare <- function(object, newdata=NULL, parameter = NULL,
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

  if (!is.null(parameter)) {
    b <- parameter
  } else {
    b <- object$coefficients
  }

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
# #' @return returns object of class "summary.glare"
# #' @export
summary.glare <- function (x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n",
                         collapse = "\n"), "\n\n", sep = "")



  coef.table <- matrix(NA, length(x$coefficients), 4L)
  dimnames(coef.table) <- list(NULL, c("Estimate",
                                       "Std. Error", "t value", "Pr(>|t|)"))
  rownames(coef.table) <- names(x$coefficients)
  coef.table[, 1] <- round(x$coefficients, digits)
  coef.table[, 2] <- round(x$coef_se, digits)
  coef.table[, 3] <- round(x$coef_z, digits)
  coef.table[, 4] <- ifelse(x$coef_p < 2e-16, "<2e-16", x$coef_p)

  cat("Coefficients:\n")
  print.default(coef.table, print.gap = 2, quote = FALSE)









  if (length(coef(x))) {
    cat("Coefficients")
    if (is.character(co <- x$contrasts))
      cat("  [contrasts: ", apply(cbind(names(co),
                                        co), 1L, paste, collapse = "="), "]")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits),
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
      x$df.residual, "Residual\n")
  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("Null Deviance:\t   ", format(signif(x$null.deviance,
                                           digits)), "\nResidual Deviance:", format(signif(x$deviance,
                                                                                           digits)), "\tAIC:", format(signif(x$aic, digits)))
  cat("\n")
  invisible(x)
}
