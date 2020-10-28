#' Anchor Generalized Linear Model main script
#'
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param A nxq matrix
#' @param xi numeric
#' @param m integer
#' @param family character
#'
#' @return
#' @export
#'
#' @examples
#' #example to be submitted
anchor_glm <- function(Y, X, A, xi, m, family){

  n <- nrow(Y)
  p <- ncol(X)

  ###############################################################
  # (to extract) Define objective likelihood
  binary_likelihood <- function(b.hat){
    return(sum(Y*(X%*%b.hat)-m*log(1+exp(X%*%b.hat))))
  }

  # (to extract) Define objective penalty
  binary_penalty <- function(b.hat){

    p.hat <- exp(X%*%b.hat)/(1+exp(X%*%b.hat))

    special.case1 <- function(Y){
      ifelse(Y==0, 0, Y*log(Y/(m*p.hat)))
    }
    special.case2 <- function(Y){
      ifelse(Y==m, 0, (m-Y)*log((m-Y)/(m-m*p.hat)))
    }

    r.D <- sign(Y/m-p.hat)*sqrt(2*(special.case1(Y)+special.case2(Y))) # deviance residuals

    P.A <- A%*%solve(t(A)%*%A)%*%t(A)
    return(t(r.D)%*%P.A%*%r.D)
  }

  ###############################################################
  # Assign objective functions depending on glm family
  likelihood_function <- switch(family,
                                binary = binary_likelihood,
                                poisson = poisson_likelihood,
                                normal = normal_likelihood
  )

  penalty_function <- switch(family,
                             binary = binary_penalty,
                             poisson = poisson_penalty,
                             normal = normal_penalty
  )

  ###############################################################
  # Construct anchor objective
  anchor_objective <- function(b.hat){
    return(1/n*(-likelihood_function(b.hat) + xi * penalty_function(b.hat)))
  }

  ###############################################################
  # Run optimization algorithm
  optimized_object <- optim(f=anchor_objective, par = runif(p), method = "L-BFGS-B")
  return(optimized_object$par)

}
