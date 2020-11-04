#' Main script for AGLM
#'
#' Script that constructs Anchor GLM objective and optimizes it
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
anchor_glm <- function(Y, X, A, xi, m = NULL, family=gaussian){

  ###############################################################
  # Initializtation
  data <- list(Y=Y, X=X, A=A)

  # Construction of model formula
  mf <- model.frame(Y ~ X-1, data)
  mm <- model.matrix(~ A-1, data)

  # Allocate glm family as in glm source code
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  # Initialize link function
  linkinv <- family$linkinv

  # Extract number of observations and parameter dimension
  n.obs <- length(Y)

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

  family$logLik <- loglike
  family$anchPen <- AnchPen

  ###############################################################
  # Construct anchor objective
  anchor_objective <- function(b.hat){
    return(1/n.obs*(-loglike(b=b.hat, Y=Y, X=X, linkinv=linkinv, m=m) + xi * AnchPen(b=b.hat, Y=Y, X=X, A=A, linkinv=linkinv, m=m)))
  }

  ###############################################################
  # Run optimization algorithm

  # Fit glm for initial parameter value
  if(family$family=="binomial" & ncol(Y)==1 & m>1){
    yy <- cbind(Y,m-Y)
  } else{
    yy <- Y
  }

  fit.glm <- glm(yy ~ X - 1, family = family$family)

  # Optimize anchor objective
  optimized_object <- stats::optim(f=anchor_objective, par = fit.glm$coefficients, method = "L-BFGS-B")

  # Construction of anchor glm class
  aglm.fit <- list(family = family,
                   m=m,
                   xi=xi,
                   mf = mf,
                   mm = mm,
                   data = data,
                   optim = optimized_object,
                   coefficients = optimized_object$par
  )
  class(aglm.fit) <- "anchorglm"

  return(aglm.fit)

}

