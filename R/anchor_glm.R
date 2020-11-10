#' Main script for AGLM
#'
#' Script that constructs Anchor GLM objective and optimizes it
#' using optim from stats.
#'
#' @param formula an object of class "formula" for the response and covariate variables.
#' @param A.formula an object of class "formula" for the anchor variables.
#' @param xi numeric
#' @param m vector
#' @param family character
#' @param type the type of residuals that should be used for the anchor penalty. The alternatives are: "deviance" (default) and "pearson". Can be abbreviated.
#'
#' @return numeric
#' @export
#' @importFrom stats glm family gaussian
anchor_glm <- function(formula, A.formula, xi, m = 1,
                       family=gaussian, type = c("deviance", "pearson")){
  # anchor_glm <- function(formula, A = NULL, xi, m = 1,
  #                        family=gaussian, type = c("deviance", "pearson")){
  ###############################################################
  # Initializtation
  type <- match.arg(type)

  # Allocate glm family as in glm source code
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  # Construction of model formula
  data <- environment(formula)
  mf <- model.frame(formula, data = data)
  Y <- model.response(mf)
  X <- model.matrix(formula, data = data)
  A <- model.matrix(A.formula, data = data)

  # Handle different form of input for binomial data
  yy <- Y # for initial parameter guess we use glm below
  if (family$family == "binomial") {
    if (dim(as.matrix(Y))[2] == 2){
      m <- yy[,1] + yy[,2]
      Y <- yy[,1]
    } else{
      yy <- cbind(Y,m-Y)
    }
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

  pearson_residuals <- function(b, Y, X, linkinv, m, family, ...){
    mu <- linkinv(X%*%b)
    V <- family$variance(mu)
    return((Y-m*mu)/sqrt(m*V))
  }

  anchor_penalty <- function(R, A, ...){
    fit <- lm(R~A)
    return(sum((fitted(fit))^2))

    # or check invertability and use
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
  if ("(Intecept)" %in% attributes(X)$dimnames[[2]]){
    glm.formula <- yy ~ X
  } else{
    glm.formula <- yy ~ X - 1
  }
  fit.glm <- glm(glm.formula, family = family$family)

  # Optimize anchor objective
  optimized_object <- stats::optim(f=anchor_objective, par = as.numeric(fit.glm$coefficients), method = "L-BFGS-B")

  # Construction of anchor glm class
  aglm.fit <- list(family = family,
                   m=m,
                   xi=xi,
                   formula = formula,
                   A.formula = A.formula,
                   data = data,
                   optim = optimized_object,
                   coefficients = optimized_object$par
  )
  class(aglm.fit) <- "anchorglm"

  return(aglm.fit)

}
