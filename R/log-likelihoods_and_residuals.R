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
binary_likelihood <- function(b, Y, X, linkinv, m, ...){

  mu <- linkinv(X%*%b)

  sum(log(choose(m,Y))+Y*log(mu)+(m-Y)*log(1-mu))
  #sum(log(choose(m,Y))+Y*(X%*%b)-m*log(1+exp(X%*%b)))
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
binary_deviance <- function(b, Y, X, linkinv, m, ...){

  mu <- linkinv(X%*%b)

  r <- Y*log(Y/(m*mu))+(m-Y)*log((m-Y)/(m-m*mu))
  p <- which(Y == 0)
  r[p] <- ((m-Y)*log((m-Y)/(m-m*mu)))[p]
  q <- which(Y == m)
  r[q] <- (Y*log(Y/(m*mu)))[q]

  sign(Y/m-mu)*sqrt(2*r)
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
poisson_likelihood <- function(b, Y, X, linkinv, ...){

  mu <- linkinv(X%*%b)

  sum(Y*log(mu)-mu)
  #sum(Y*log(mu)-mu-log(factorial(Y)))

  #sum(Y*(X%*%b)-exp(X%*%b))
  #sum(Y*(X%*%b)-exp(X%*%b)-log(factorial(Y)))
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
poisson_deviance <- function(b, Y, X, linkinv, ...){

  mu <- linkinv(X%*%b)

  r <- mu
  p <- which(Y > 0)
  r[p] <- (Y * log(Y/mu) - (Y - mu))[p]

  sign(Y-mu)*sqrt(2 * r)
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
normal_likelihood <- function(b, Y, X, linkinv, ...){

  #n.obs <- nrow(Y)

  mu <- linkinv(X%*%b)
  #ss <- sum((Y-mu)^2)/n.obs

  #-n.obs/2*log(2*pi*ss)-1/(2*ss)*sum((Y-mu)^2)

  #-1/2*sum((Y-mu)^2)
  -sum((Y-mu)^2)
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
normal_deviance <- function(b, Y, X, linkinv, ...){

  mu <- linkinv(X%*%b)

  Y - mu
}
