#' Binary log likelihood
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param m integer
#'
#' @return log likelihood objective
#' @export
binary_likelihood <- function(b, Y, X, m){
  return(sum(Y*(X%*%b)-m*log(1+exp(X%*%b))))
}

#' Binary anchor penalty
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param m integer
#' @param P.A orthogonal projection matrix of A
#'
#' @return anchor penalty objective
#' @export
binary_penalty <- function(b, Y, X, m, P.A){

  p.hat <- exp(X%*%b)/(1+exp(X%*%b))

  special.case1 <- function(Y){
    ifelse(Y==0, 0, Y*log(Y/(m*p.hat)))
  }
  special.case2 <- function(Y){
    ifelse(Y==m, 0, (m-Y)*log((m-Y)/(m-m*p.hat)))
  }

  r.D <- sign(Y/m-p.hat)*sqrt(2*(special.case1(Y)+special.case2(Y))) # deviance residuals

  return(t(r.D)%*%P.A%*%r.D)
}

##################################################################################

#' Poisson likelihood
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param m integer
#'
#' @return log likelihood objective
#' @export
poisson_likelihood <- function(b, Y, X, m){
  return(sum(Y*(X%*%b)-exp(X%*%b)))
}

#' Poisson anchor penalty
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param m integer
#' @param P.A orthogonal projection matrix of A
#'
#' @return anchor penalty objective
#' @export
poisson_penalty <- function(b, Y, X, m, P.A){

  mu <- exp(X%*%b) # inverse of logit link

  special.case <- function(Y,mu){
    ifelse(Y==0|mu==0, mu, Y*log(Y/mu)-(Y-mu))
  }

  r.D <- sign(Y-mu)*sqrt(2*special.case(Y,mu))
  return(t(r.D)%*%P.A%*%r.D)
}

##################################################################################

#' Normal likelihood
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param m integer
#'
#' @return log likelihood objective
#' @export
normal_likelihood <- function(b, Y, X, m){

  n.obs <- nrow(Y)

  mu <- X%*%b
  ss <- sum((Y-mu)^2)/n.obs

  return(-n.obs/2*log(2*pi*ss)-1/(2*ss)*sum((Y-mu)^2))
}

#' Normal anchor penalty
#'
#' @param b parameter vector
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param m integer
#' @param P.A orthogonal projection matrix of A
#'
#' @return anchor penalty objective
#' @export
normal_penalty <- function(b, Y, X, m, P.A){

  mu <- X%*%b

  r.D <- Y - mu
  return(t(r.D)%*%P.A%*%r.D)
}
