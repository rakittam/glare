#' Main script for AGLM
#'
#' SCript that constructs Anchor GLM objective and optimizes it
#' using optim from stats.
#'
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param A nxq matrix
#' @param xi numeric
#' @param m integer
#' @param family character
#'
#' @return numeric
#' @export
#' @examples
#' Y <- c(0,1,0)
#' X <- matrix(c(1,2,3, 11,12,13), nrow = 3, ncol = 2)
#' A <- matrix(c(1,-1,-1, 1,-1,1), nrow = 3, ncol = 2)
#' anchor_glm(Y, X, A, 2, 1, "binomial")
#' @importFrom stats glm family gaussian
anchor_glm <- function(Y, X, A, xi, m = 0, family=gaussian){

  ###############################################################
  # Initializtation

  # Allocate glm family as in glm source code
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  # Extract number of observations and parameter dimension
  n <- length(Y)
  p <- ncol(X)

  # Calculate orthogonal projection onto column space of A
  P.A <- A%*%solve(t(A)%*%A)%*%t(A)

  ###############################################################
  # Define needed functions

  # Initialize weights, as needed for input in family
  weights <- rep.int(1, n)

  # Extract aic and deviance residual functions from family
  aic <- family$aic
  dev.resids <- family$dev.resids

  # Construct loglikelihood and anchor penalty
  loglike <- function(b, Y, X, m, weights, p){

    mu <- X%*%b
    dev <- sum(dev.resids(Y, mu, weights))

    return(1/2*aic(Y, m, mu, weights, dev) - p)
  }

  AnchPen <- function(b, Y, X, P.A, weights){

    mu <- X%*%b
    r.D <- dev.resids(Y, mu, weights)

    return(t(r.D)%*%P.A%*%r.D)
  }

  ###############################################################
  # Construct anchor objective
  anchor_objective <- function(b.hat){
    return(1/n*(-loglike(b.hat, Y, X, m, weights, p) + xi * AnchPen(b.hat, Y, X, P.A, weights)))
  }

  ###############################################################
  # Run optimization algorithm

  # Fit glm for initial parameter guess
  #fit.glm <- glm(Y, X, family)

  optimized_object <- stats::optim(f=anchor_objective, par = stats::runif(p), method = "L-BFGS-B")

  aglm.fit <- list(par = optimized_object$par,
                   value = optimized_object$value,
                   family = family)
  class(aglm.fit) <- "anchorglm"

  return(aglm.fit)

}
