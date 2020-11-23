#' Data extraction
#'
#' Returns data from object input.
#'
#' @param object an object of class \code{"\link{glare}"}.
#' @param newdata an optional data frame, list or environment containing the
#'  variables in the model. If not used, the variables are taken from the
#'  objective.
#' @param parameter optional parameter input. Default takes the parameter
#'  from the glare objective.
#' @param ... further arguments passed to or from other methods.
#'
#' @return list of data extracted from input.
#' @export
#' @importFrom stats model.response model.matrix model.frame
extract_data <- function(object, newdata = NULL, parameter = NULL, ...) {

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

  list(Y = Y, X = X, b = b, m = m)
}

#' Extract Log-Likelihood
#'
#' Returns log-likelihood of object with class `"glare"`, using initial
#'  training or new data.
#'
#' @param object an object of class \code{"\link{glare}"}.
#' @param newdata an optional data frame, list or environment containing the
#'  variables in the model. If not used, the variables are taken from the
#'  objective.
#' @param parameter optional parameter input. Default takes the parameter
#'  from the glare objective.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Returns the log-likelihood of an object of class `"glare"`.
#' @export
logLik.glare <- function(object, newdata = NULL, parameter = NULL, ...) {

  data <- extract_data(object = object,
                       newdata = newdata,
                       parameter = parameter)

  object$logLik(b = data$b,
                Y = data$Y,
                X = data$X,
                linkinv = object$family$linkinv,
                m = data$m)
}

#' Extract Model Coefficients
#'
#' Returns model coefficients of object with class `"glare"`.
#'
#' @param object an object of class \code{"\link{glare}"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Coefficients extracted from the model object of class `"glare"`.
#' @export
coef.glare <- function(object, ...) {
  object$coefficients
}

#' Model Predictions
#'
#' Returns predictions on scale of linear predictor or response of `"glare"`
#'  object using train or new data.
#'
#' @param object an object of class \code{"\link{glare}"}.
#' @param newdata an optional data frame, list or environment containing the
#'  variables in the model. If not used, the variables are taken from the
#'  objective.
#' @param parameter optional parameter input. Default takes the parameter
#'  from the glare objective.
#' @param type the scale of prediction. The default is on the scale of
#'  the linear predictors; the alternative `"response"` is on the scale of the
#'  response variable. Thus for a default binomial model the default predictions
#'  are of log-odds (probabilities on logit scale) and `type = "response"` gives
#'  the predicted probabilities.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The returned predictions depend on the choosen type of scale (see
#'  above).
#' @export
predict.glare <- function(object, newdata = NULL, parameter = NULL,
                              type = c("link", "response"), ...) {

  type <- match.arg(type)

  data <- extract_data(object = object,
                       newdata = newdata,
                       parameter = parameter)

  switch(type,
         "link"= data$X %*% data$b,
         "response"= object$family$linkinv(data$X %*% data$b))
}

#' Extract Model Residuals
#'
#' Returns deviance or pearson residuals of a `"glare"` object using train or
#'  new data.
#'
#' @param object an object of class \code{"\link{glare}"}.
#' @param newdata an optional data frame, list or environment containing the
#'  variables in the model. If not used, the variables are taken from the
#'  objective.
#' @param parameter optional parameter input. Default takes the parameter
#'  from the glare objective.
#' @param type the type of residuals which should be returned. The alternatives
#'  are: "deviance" (default) and "pearson". Can be abbreviated.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Residuals extracted from the model object of class `"glare"`.
#' @export
residuals.glare <- function(object, newdata=NULL, parameter = NULL,
                                type = c("deviance", "pearson"), ...) {

  type <- match.arg(type)

  data <- extract_data(object = object,
                       newdata = newdata,
                       parameter = parameter)

  res_function <- switch(type,
                         "deviance" = object$devianceRes,
                         "pearson" = object$pearsonRes)

  res_function(b = data$b, Y = data$Y, X = data$X,
               linkinv = object$family$linkinv,
               m = data$m, family = object$family)
}


#' Object Summary
#'
#' Returns summary of an object of class `"glare"`.
#'
#' @param object an object of class \code{"\link{glare}"}.
#' @param digits integer for number of digits used for the rounded result.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Returns object of class `"summary.glare"`.
#' @export
summary.glare <- function (object, digits = max(3L, getOption("digits") - 3L), ...) {

  # Deviance Residuals
  res.table <- matrix(NA, 1, 6L)
  dimnames(res.table) <- list(NULL, c("Min.", "1st Qu.", "Median",
                                      "Mean", "3rd Qu.", "Max."))
  rownames(res.table) <- ""
  res.table[,1] <- round(min(object$r_D), digits)
  res.table[1,2] <- as.numeric(round(quantile(object$r_D)[2], digits))
  res.table[1,3] <- round(median(object$r_D), digits)
  res.table[1,4] <- round(mean(object$r_D), digits)
  res.table[1,5] <- as.numeric(round(quantile(object$r_D)[4], digits))
  res.table[1,6] <- round(max(object$r_D), digits)

  # Coefficient estimates
  coef.table <- matrix(NA, length(object$coefficients), 4L)
  dimnames(coef.table) <- list(NULL, c("Estimate",
                                       "Std. Error", "t value", "Pr(>|t|)"))
  rownames(coef.table) <- names(object$coefficients)
  coef.table[, 1] <- round(object$coefficients, digits)
  coef.table[, 2] <- round(object$coef_se, digits)
  coef.table[, 3] <- round(object$coef_z, digits)
  coef.table[, 4] <- ifelse(object$coef_p < 2e-16, "<2e-16", object$coef_p)

  # Output summary class
  ans <- list(fit = object,
              digits = digits,
              res.table = res.table,
              coef.table = coef.table)
  class(ans) <- "summary.glare"
  return(ans)
}

#' Print of summary method for `"glare"` object.
#'
#' Returns summary of a glare object.
#'
#' @param object an object of class \code{"\link{glare}"}.
#' @param digits integer for number of digits used for the rounded result.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Prints summary of a `"glare"` object.
#' @export
print.summary.glare <- function (object, digits = max(3L, getOption("digits") - 3L), ...) {

  summary_glare <- summary.glare(object, digits = digits)

  # Function call
  cat("\nCall:  ", paste(deparse(object$call), sep = "\n",
                         collapse = "\n"), "\n\n", sep = "")

  # Deviance Residuals
  cat("Deviance Residuals:\n")
  print.default(summary_glare$res.table, print.gap = 2, quote = FALSE)
  cat("\n")

  # Coefficient estimates
  cat("Coefficients:\n")
  print.default(summary_glare$coef.table, print.gap = 2, quote = FALSE)
  cat("\n")

  invisible(object)
}



#
# summary.glare <- function(object, ...) {
#   print.summary.glare(object)
# }
#
# summary <- function(object, ...) {
#   UseMethod("summary")
# }
