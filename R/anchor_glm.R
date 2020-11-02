#' Main script for AGLM
#'
#' SCript that constructs Anchor GLM objective and optimizes it
#' using optim from stats.
#'
#' @param Y n-dimensional vector
#' @param X nxp matrix
#' @param A nxq matrix
#' @param xi numeric
#' @param m vector
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
anchor_glm <- function(Y, X, A, xi, m = 1, family=gaussian){

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

  #############################################################
  # Assign objective functions depending on glm family
  loglike <- switch(family$family,
                    "binomial" = binary_likelihood,
                    "poisson" = poisson_likelihood,
                    "gaussian" = normal_likelihood
  )

  AnchPen <- switch(family$family,
                    "binomial" = binary_penalty,
                    "poisson" = poisson_penalty,
                    "gaussian" = normal_penalty
  )

  ###############################################################
  # Construct anchor objective
    anchor_objective <- function(b.hat){
    return(1/n*(-loglike(b.hat, Y, X, m) + xi * AnchPen(b.hat, Y, X, m, P.A)))
  }

  ###############################################################
  # Run optimization algorithm

  # Fit glm for initial parameter guess

  yy <- Y
  if(family$family=="binomial"){
    yy <- as.factor(Y)
  }
  fit.glm <- glm(yy ~ X - 1, family = family$family)

  #optimized_object <- stats::optim(f=anchor_objective, par = stats::runif(p), method = "L-BFGS-B")
  optimized_object <- stats::optim(f=anchor_objective, par = fit.glm$coefficients, method = "L-BFGS-B")

  aglm.fit <- list(par = optimized_object$par,
                   value = optimized_object$value,
                   family = family)
  class(aglm.fit) <- "anchorglm"

  return(aglm.fit)

}
