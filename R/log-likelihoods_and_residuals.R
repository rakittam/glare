#' Binary log likelihood
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param linkinv inversion of link function
#' @param m integer
#' @param indiv logical input if log-likelihood of each observation or the sum
#'  of the log-likelihoods should be returned.
#' @param ... other
#'
#' @return log likelihood objective
#' @export
#' @importFrom stats dbinom dnorm dpois
binary_likelihood <- function(b, Y, X, linkinv, m, indiv = FALSE, ...) {

  mu <- linkinv(X %*% b)
  logLik_indiv <- dbinom(Y, m, mu, log = TRUE)

  if (indiv) {
    return(logLik_indiv)
  } else {
    return(sum(logLik_indiv))
  }
}

#' Binary deviance residuals
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param linkinv inversion of link function
#' @param m integer
#' @param ... other
#'
#' @return deviance residuals
#' @export
binary_deviance <- function(b, Y, X, linkinv, m, ...) {

  mu <- linkinv(X %*% b)

  r <- Y * log(Y / (m * mu)) + (m - Y) * log((m - Y) / (m - m * mu))
  p <- which(Y == 0)
  r[p] <- ((m - Y) * log((m - Y) / (m - m * mu)))[p]
  q <- which(Y == m)
  r[q] <- (Y * log(Y / (m * mu)))[q]

  sign(Y / m - mu) * sqrt(2 * r)
}

# -------------------------------------------------------------------
#' Poisson likelihood
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param linkinv inversion of link function
#' @param indiv logical input if log-likelihood of each observation or the sum
#'  of the log-likelihoods should be returned.
#' @param ... other
#'
#' @return log likelihood objective
#' @export
poisson_likelihood <- function(b, Y, X, linkinv, indiv = FALSE, ...) {

  mu <- linkinv(X %*% b)
  logLik_indiv <- dpois(Y, mu, log = TRUE)

  if (indiv) {
    return(logLik_indiv)
  } else {
    return(sum(logLik_indiv))
  }
}

#' Poisson deviance residuals
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param linkinv inversion of link function
#' @param ... other
#'
#' @return deviance residuals
#' @export
poisson_deviance <- function(b, Y, X, linkinv, ...) {

  mu <- linkinv(X %*% b)

  r <- mu
  p <- which(Y > 0)
  r[p] <- (Y * log(Y / mu) - (Y - mu))[p]

  sign(Y - mu) * sqrt(2 * r)
}

# ----------------------------------------------------------------------
#' Normal likelihood
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param linkinv inversion of link function
#' @param indiv logical input if log-likelihood of each observation or the sum
#'  of the log-likelihoods should be returned.
#' @param ... other
#'
#' @return log likelihood objective
#' @export
normal_likelihood <- function(b, Y, X, linkinv, indiv = FALSE, ...) {

  mu <- linkinv(X %*% b)
  n <- length(Y)
  s <- sqrt(1/n * sum((Y-mu)^2))
  logLik_indiv <- dnorm(Y, mean = mu, sd = s, log = TRUE)

  if (indiv) {
    return(logLik_indiv)
  } else {
    return(sum(logLik_indiv))
  }
}

#' Normal deviance residuals
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param linkinv inversion of link function
#' @param ... other
#'
#' @return deviance residuals
#' @export
normal_deviance <- function(b, Y, X, linkinv, ...) {

  mu <- linkinv(X %*% b)
  n <- length(Y)
  s <- sqrt(1/n * sum((Y-mu)^2))
  1 / s * (Y - mu)
}

#' Normal classical residuals
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param linkinv inversion of link function
#' @param ... other
#'
#' @return deviance residuals
#' @export
normal_classical <- function(b, Y, X, linkinv, ...) {

  mu <- linkinv(X %*% b)
  Y - mu
}
