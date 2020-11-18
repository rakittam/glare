#' Binary log likelihood
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param linkinv inversion of link function
#' @param m integer
#' @param ... other
#'
#' @return log likelihood objective
#' @export
#' @importFrom stats dbinom dnorm dpois
binary_likelihood <- function(b, Y, X, linkinv, m, ...) {

  mu <- linkinv(X %*% b)

  #sum(log(choose(m, Y)) + Y * log(mu) + (m - Y) * log(1 - mu))
  sum(dbinom(Y, m, mu, log = TRUE))
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
#' @param ... other
#'
#' @return log likelihood objective
#' @export
poisson_likelihood <- function(b, Y, X, linkinv, ...) {

  mu <- linkinv(X %*% b)

  #sum(Y * log(mu) - mu)
  #sum(Y * log(mu) - mu - log(factorial(Y)))
  sum(dpois(Y, mu, log = TRUE))
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
#' @param ... other
#'
#' @return log likelihood objective
#' @export
normal_likelihood <- function(b, Y, X, linkinv, ...) {

  mu <- linkinv(X %*% b)

  #-sum((Y - mu)^2)

  n <- length(Y)
  s <- sqrt(1/n * sum((Y-mu)^2))
  sum(dnorm(Y, mean = mu, sd = s, log = TRUE))
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

  Y - mu
}
