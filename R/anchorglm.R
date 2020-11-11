#' Fitting Anchored Generalized Linear Models
#'
#' `anchorglm` is used to fit anchored generalized linear models, specified by giving an anchor variable, a symbolic description of the linear predictor and a description of the error distribution.
#'
#' @param formula an object of class \code{"\link{formula}"} for the response and covariate variables.
#' @param A.formula an object of class \code{"\link{formula}"} for the anchor variables.
#' @param data an optional data frame, list or environment containing the variables in the model. If not found in data, the variables are taken from `environment(formula)`.
#' @param xi a numeric value for the hyperparameter xi.
#' @param m number of independent bernoulli trials. See the vignette with the binomial example for help.
#' @param family a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. See \code{\link{family}} for details. The families supported up to now are `gaussian`, `binomial` and `poisson`.
#' @param type the type of residuals that should be used for the anchor penalty. The supported alternatives are `"deviance"` (default) and `"pearson"`. Can be abbreviated.
#'
#' @details The response-covariate formula must have the form `response ~ covariates`, where `response` is the response vector and `covariates` is a series of terms which specifies a linear predictor for the response.\cr\cr
#' For binomial families the response can also be specified as a two-column matrix with the columns giving the numbers of successes and failures.\cr\cr
#' The anchor formula must have the form `~ anchors`, where `anchors` contains all given anchor variables.
#'
#' @return `anchorglm` returns an object of class `"anchorglm"`.\cr\cr
#' The function \code{\link{summary}} (i.e., \code{\link{summary.anchorglm}}) can be used to obtain or print a summary of the results.\cr\cr
#' The generic accessor functions `logLik`, `coef`, `predict` and `residuals` can be used to extract various useful features of the value returned by `anchorglm`.\cr\cr
#' An object of class `"anchorglm"` is a list containing at least the following components:
#' \item{formula}{the response and covariate \code{"\link{formula}"} used.}
#' \item{A.formula}{the anchor \code{"\link{formula}"} used.}
#' \item{data}{the `data` argument.}
#' \item{xi}{the hyperparameter used.}
#' \item{m}{the used bernoulli trials coefficients. If not used set to `1`.}
#' \item{family}{the \code{"\link{formula}"} object used.}
#' \item{logLik}{the log-likelihood function used corresponding to the used family.}
#' \item{devianceRes}{the deviance residuals function used corresponding to the used family.}
#' \item{pearsonRes}{the pearson residuals function used corresponding to the used family.}
#' \item{optim}{the \code{\link{optim}} output as a list.}
#' \item{coefficients}{the optimal coefficients derived by optim.}
#' @export
#' @importFrom stats glm family gaussian model.response model.matrix model.frame lm fitted optim
anchorglm <- function(formula, A.formula, data, xi, m = 1,
                       family=gaussian, type = c("deviance", "pearson")){
  ###############################################################
  # Initializtation
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
  A <- model.matrix(A.formula, data = data)

  # Handle different form of input for binomial data
  yy <- Y # for initial parameter guess for optim we use glm, see below
  if (family$family == "binomial") {
    if (dim(as.matrix(Y))[2] == 2){
      m <- yy[,1] + yy[,2]
      Y <- yy[,1]
    } else{
      yy <- cbind(Y,m-Y)
    }
  }

  # Call link function
  linkinv <- family$linkinv

  # Extract number of observations
  n.obs <- length(Y)

  #############################################################
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

  pearson_residuals <- function(b, Y, X, linkinv, m, family, ...){
    mu <- linkinv(X%*%b)
    V <- family$variance(mu)
    return((Y-m*mu)/sqrt(m*V))
  }

  anchor_penalty <- function(R, A, ...){
    fit.temp <- lm(R~A)
    return(sum((fitted(fit.temp))^2))
  }

  if(type == "deviance"){
    anchRes <- deviance_residuals
  }
  if(type == "pearson"){
    anchRes <- pearson_residuals
  }

  ###############################################################
  # Construct anchor objective
  anchor_objective <- function(b.hat){
    return(1/n.obs*(-log_likelihood(b=b.hat, Y=Y, X=X, linkinv=linkinv, m=m) +
                      xi * anchor_penalty(R=anchRes(b=b.hat, Y=Y, X=X, linkinv=linkinv, m=m, family=family), A=A)))
  }

  ###############################################################
  # Run optimization algorithm

  # Fit glm for initial parameter value
  if ("(Intecept)" %in% attributes(X)$dimnames[[2]]){
    glm.formula <- yy ~ X
  } else{
    glm.formula <- yy ~ X - 1
  }
  fit.glm <- glm(glm.formula, family = family$family)

  # Optimize anchor objective
  optimized_object <- optim(fn=anchor_objective, par = as.numeric(fit.glm$coefficients), method = "L-BFGS-B", hessian = TRUE)

  # Construction of anchor glm class
  aglm.fit <- list(formula = formula,
                   A.formula = A.formula,
                   data = data,
                   xi=xi,
                   m=m,
                   family = family,
                   logLik = log_likelihood,
                   devianceRes = deviance_residuals,
                   pearsonRes = pearson_residuals,
                   optim = optimized_object,
                   coefficients = optimized_object$par
  )
  class(aglm.fit) <- "anchorglm"

  return(aglm.fit)
}
