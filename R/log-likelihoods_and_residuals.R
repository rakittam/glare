##################################################################################

#' Binary log likelihood
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param m integer
#' @param ... other
#'
#' @return log likelihood objective
#' @export
binary_likelihood <- function(b, Y, X, m, ...){
  return(sum(Y*(X%*%b)-m*log(1+exp(X%*%b))))
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

  r.D <- sign(Y/m-mu)*sqrt(2*r) # deviance residuals

  return(r.D)
}

##################################################################################

#' Poisson likelihood
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param ... other
#'
#' @return log likelihood objective
#' @export
poisson_likelihood <- function(b, Y, X, ...){
  return(sum(Y*(X%*%b)-exp(X%*%b)))
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
  r.D <- sign(Y-mu)*sqrt(2 * r)

  return(r.D)
}

##################################################################################

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

  n.obs <- nrow(Y)

  mu <- linkinv(X%*%b)
  ss <- sum((Y-mu)^2)/n.obs

  return(-n.obs/2*log(2*pi*ss)-1/(2*ss)*sum((Y-mu)^2))
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
  r.D <- Y - mu

  return(r.D)
}
