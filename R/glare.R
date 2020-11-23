#' Fitting Generalized Linear Anchor Regression
#'
#' `glare` is used to fit Generalized Linear Anchor REgression (GLARE), specified by
#'  giving an anchor variable, a symbolic description of the linear predictor
#'  and a description of the error distribution.
#'
#' @param formula an object of class \code{"\link{formula}"} for the response
#'  and covariate variables.
#' @param A_formula an object of class \code{"\link{formula}"} for the anchor
#'  variables.
#' @param data an optional data frame, list or environment containing the
#'  variables in the model. If not found in data, the variables are taken from
#'   `environment(formula)`.
#' @param xi a numeric value for the hyperparameter xi.
#' @param family a description of the error distribution and link function to
#'  be used in the model. This can be a character string naming a family
#'  function, a family function or the result of a call to a family function.
#'  See \code{\link{family}} for details. The families supported up to now are
#'  `gaussian`, `binomial` and `poisson`.
#' @param type the type of residuals that should be used for the anchor
#'  penalty. The supported alternatives are `"deviance"` (default) and
#'  `"pearson"`. Can be abbreviated.
#'
#' @details The response-covariate formula must have the form
#'  `response ~ covariates`, where `response` is the response vector and
#'  `covariates` is a series of terms which specifies a linear predictor for
#'  the response.\cr\cr
#'  For binomial families the response can also be specified as a two-column
#'  matrix with the columns giving the numbers of successes and failures. If the
#'  number of independent bernoulli trials is provided trough the data
#'  argument, the variable has be called `m`. See the vignette with the
#'  binomial example for help.\cr\cr
#'  The anchor formula must have the form `~ anchors`, where `anchors` contains
#'  all given anchor variables.
#'
#' @return `glare` returns an object of class `"glare"`.\cr\cr
#'  The function \code{\link{summary}} (i.e., \code{\link{summary.glare}})
#'  can be used to obtain or print a summary of the results.\cr\cr
#'  The generic accessor functions `logLik`, `coef`, `predict` and `residuals`
#'  can be used to extract various useful features of the value returned by
#'  `glare`.\cr\cr
#'  An object of class `"glare"` is a list containing at least the
#'  following components:
#'  \item{formula}{the response and covariate \code{"\link{formula}"} used.}
#'  \item{A_formula}{the anchor \code{"\link{formula}"} used.}
#'  \item{data}{the `data` argument as data.frame.}
#'  \item{xi}{the hyperparameter used.}
#'  \item{family}{the \code{"\link{formula}"} object used.}
#'  \item{logLik}{the log-likelihood function used corresponding to the used
#'  family.}
#'  \item{devianceRes}{the deviance residuals function used corresponding to
#'  the used family.}
#'  \item{pearsonRes}{the pearson residuals function used corresponding to the
#'  used family.}
#'  \item{optim}{the \code{\link{optim}} output as a list.}
#'  \item{coefficients}{the optimal coefficients derived by optim.}
#' @export
#' @importFrom stats glm family gaussian model.response model.matrix model.frame
#'  lm fitted optim pnorm median quantile
glare <- function(formula, A_formula, data, xi,
                       family = gaussian, type = c("deviance", "pearson")) {
  # Initializtation -----------------------------------------------
  cal <- match.call()
  type <- match.arg(type)

  # Allocate glm family as in glm source code
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  # Construction of model formula
  if (missing(data))
    data <- environment(formula)
  mf <- model.frame(formula, data = data)
  Y <- model.response(mf)
  X <- model.matrix(formula, data = data)
  A <- model.matrix(A_formula, data = data)

  Xnames <- dimnames(X)[[2L]]
  Yname <- dimnames(mf)[[2]][1]

  # Handle different form of input for binomial data
  yy <- Y # for initial parameter guess for optim we use glm, see below

  if (family$family == "binomial") {
    if (dim(as.matrix(Y))[2] == 2) {
      m <- yy[, 1] + yy[,2]
      Y <- yy[, 1]
    } else if ("m" %in% names(data)) {
      m <- data$m
      yy <- cbind(Y, m - Y)
    }
    else {
      m <- 1
      yy <- cbind(Y, m - Y)
    }
  } else {
    m <- 1
  }

  # Call link function
  linkinv <- family$linkinv

  # Extract number of observations
  n.obs <- length(Y)

  # Function assignment -----------------------------------------------
  # Assign needed functions depending on given glm family
  log_likelihood <- switch(family$family,
                           "binomial" = binary_likelihood,
                           "poisson" = poisson_likelihood,
                           "gaussian" = normal_likelihood
  )

  deviance_residuals <- switch(family$family,
                               "binomial" = binary_deviance,
                               "poisson" = poisson_deviance,
                               "gaussian" = normal_deviance
  )

  pearson_residuals <- function(b, Y, X, linkinv, m, family, ...) {
    mu <- linkinv(X %*% b)
    V <- family$variance(mu)
    (Y - m * mu)/sqrt(m * V)
  }

  anchor_penalty <- function(R, A, ...) {

    #P.A <- A %*% solve(t(A) %*% A) %*% t(A)
    #return(t(R)%*%P.A%*%R)

    fit_temp <- lm(R ~ A)
    sum((fitted(fit_temp))^2)
  }

  if (type == "deviance") {
    anch_res <- deviance_residuals
  }
  if (type == "pearson") {
    anch_res <- pearson_residuals
  }

  # Construction of anchor objective ------------------------------------------
  anchor_objective <- function(b_hat) {
    1/n.obs * (-log_likelihood(b = b_hat,
                               Y = Y,
                               X = X,
                               linkinv = linkinv,
                               m = m) +
               xi * anchor_penalty(R = anch_res(b = b_hat,
                                                Y = Y,
                                                X = X,
                                                linkinv = linkinv,
                                                m = m,
                                                family = family),
                                   A = A))
  }

  # Run optimization algorithm -----------------------------------------------

  # Fit glm for initial parameter value
  if ("(Intecept)" %in% attributes(X)$dimnames[[ 2 ]]) {
    glm_formula <- yy ~ X
  } else {
    glm_formula <- yy ~ X - 1
  }
  glm_fit <- glm(glm_formula, family = family$family)

  # Optimize anchor objective
  optimized_object <- optim(fn = anchor_objective,
                            par = as.numeric(glm_fit$coefficients),
                            method = "L-BFGS-B",
                            hessian = TRUE)

  # Coefficients
  coefficients <- optimized_object$par
  names(coefficients) <- Xnames
  hessian <- optimized_object$hessian
  coef_se <- diag(sqrt(solve(hessian)*1/n.obs))
  coef_z <- coefficients/coef_se
  coef_p <- pnorm(-abs(coef_z))  * 2

  # Deviance Residuals:
  r_D <- deviance_residuals(coefficients, Y, X, linkinv, m)

  # Construction of anchor glm class
  glare_fit <- list(call = cal,
                    formula = formula,
                    A_formula = A_formula,
                    data = data,
                    xi = xi,
                    m = m,
                    family = family,
                    logLik = log_likelihood,
                    devianceRes = deviance_residuals,
                    pearsonRes = pearson_residuals,
                    optim = optimized_object,
                    coefficients = coefficients,
                    coef_se = coef_se,
                    coef_z = coef_z,
                    coef_p = coef_p,
                    r_D = r_D
  )
  class(glare_fit) <- "glare"

  glare_fit
}
