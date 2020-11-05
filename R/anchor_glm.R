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
#' @param type the type of residuals that should be used for the anchor penalty. The alternatives are: "deviance" (default) and "pearson". Can be abbreviated.
#'
#' @return numeric
#' @export
#' @importFrom stats glm family gaussian
anchor_glm <- function(Y, X, A, xi, m = 1,
                       family=gaussian, type = c("deviance", "pearson")){

  type <- match.arg(type)

  # HERE einbauen, if family binomial -> check for form of input and generalize it for this script
  # use PDF on desktop

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

  # pearson_residuals <- switch(family$family,
  #                       "binomial" = binary_pearson,
  #                       "poisson" = poisson_pearson,
  #                       "gaussian" = normal_pearson
  # )

  pearson_residuals <- function(b, Y, X, linkinv, m, family, ...){
    mu <- linkinv(X%*%b)
    V <- family$variance(mu)
    return((Y-m*mu)/sqrt(m*V))
  }

  anchor_penalty <- function(R, A, ...){
    fit <- lm(R~A)
    return(sum((fitted(fit))^2))

    #P.A <- A%*%solve(t(A)%*%A)%*%t(A)
    #return(t(r.D)%*%P.A%*%r.D)
  }

  family$logLik <- log_likelihood
  family$devianceRes <- deviance_residuals
  family$pearsonRes <- pearson_residuals
  family$anchPen <- anchor_penalty

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

