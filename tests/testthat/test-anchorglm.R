test_that("Testing construction of binomial anchor objective and optimization", {

  set.seed(1)

  n <- 1000 # number of samples from unpertubed and pertubed distribution

  # Initialize data

  # Anchor coefficients
  g1 <- 0.5
  g2 <- -0.2
  g3 <- -0.4
  g4 <- -2

  m <- sample(1:5, size = n, replace = TRUE) # number of trials for binary distribution

  # Initialize training data
  A <- matrix(sample(c(-1, 1), size = n * 2, replace = TRUE), nrow = n, ncol = 2)

  epsH <- matrix(stats::rnorm(n = n * 1, mean = 0, sd = 1), nrow = n, ncol = 1)
  H <- epsH

  epsX <- matrix(stats::rnorm(n = n * 2, mean = 0, sd = 1), nrow = n, ncol = 2)
  X <- matrix(nrow = n, ncol = 2)
  X[, 1] <- g1*A[, 1]+g2*A[, 2]+H+epsX[, 1]
  X[, 2] <- g1*A[, 1]+g3*A[, 2]+H+epsX[, 2]

  Y <- matrix(stats::rbinom(n=n, size=m, stats::plogis(3*X[,1]+3*X[,2]+H+g4*A[,1])), nrow = n, ncol = 1)

  # Set up test bench
  P.A <- A%*%solve(t(A)%*%A)%*%t(A)

  AGLM <- function(xi){

    # Step 1. Define the objective loss
    loss <- function(b.hat){
      return(-sum(log(choose(m,Y))+Y*(X%*%b.hat)-m*log(1+exp(X%*%b.hat))))
    }

    # Step 2. Define anchor penalty
    anchor_penalty <- function(b.hat){
      p.hat <- exp(X%*%b.hat)/(1+exp(X%*%b.hat)) # inverse of logit link

      special.case1 <- function(Y){
        ifelse(Y==0, 0, Y*log(Y/(m*p.hat)))
      }
      special.case2 <- function(Y){
        ifelse(Y==m, 0, (m-Y)*log((m-Y)/(m-m*p.hat)))
      }

      r.D <- sign(Y/m-p.hat)*sqrt(2*(special.case1(Y)+special.case2(Y))) # deviance residuals

      fit.temp <- lm(r.D~A)
      return(sum((fitted(fit.temp))^2))
      #return(t(r.D)%*%P.A%*%r.D)
    }

    # Step 3. Contruct objective by 1. and 2.
    objective <- function(b.hat){
      return(1/n*(loss(b.hat) + xi * anchor_penalty(b.hat)))
    }

    # Set start value for optimization
    ans2 <- optim(f=objective, par = runif(2), method = "L-BFGS-B")
    return(ans2$par)

    #ans2 <- optim(f=objective, par = runif(ncol(X)), method = "L-BFGS-B")
  }

  xi=0
  YY <- cbind(Y, m-Y)
  fit.glm <- glm(formula = YY~X-1, family = binomial)
  fit.aglm <- anchorglm(formula = YY~X-1, A_formula = ~A-1, xi=xi, m=m, family=binomial, type="deviance")
  logLik(fit.glm)
  logLik(fit.aglm)

  xi=2
  AGLM(xi=xi)
  as.numeric(anchorglm(formula = Y~X-1, A_formula = ~A-1, xi=xi, m=m, family=binomial, type="deviance")$optim$par)
  as.numeric(anchorglm(formula = Y~X-1, A_formula = ~A-1, xi=xi, m=m, family=binomial, type="pearson")$optim$par)

  # fit.aglm <- anchorglm(formula = Y~X-1, A_formula = ~A-1, xi=xi, m=m, family=binomial, type="deviance")
  #
  # fit.aglm2 <- anchorglm(formula = Y~X-1, A_formula = ~A-1, xi=xi, m=m, family=binomial, type="deviance")
  #
  #
  # # Or with (success, fails)
  # YY <- cbind(Y, m-Y)
  # as.numeric(anchorglm(formula = YY~X-1, A_formula = ~A-1, xi=xi, family=binomial, type="pearson")$optim$par)
  # fit.glm <- glm(formula = YY~X-1, family = binomial)
  # #fit.aglm3 <-
  # logLik(fit.aglm)
  # logLik(fit.glm)
  # logLik(fit.aglm2)

  # Compare results
  expect_equal(AGLM(2),
               as.numeric(anchorglm(formula = Y~X-1, A_formula = ~A-1, xi=xi, m=m, family=binomial, type="deviance")$optim$par), tolerance = 0.00001)
})

#####################################################################################
test_that("Testing construction of poisson anchor objective and optimization", {

  set.seed(1)

  n <- 1000 # number of samples from unpertubed and pertubed distribution

  # Initialize data

  # Anchor coefficients
  g1 <- 0.5

  m <- sample(1:5, size = n, replace = TRUE) # number of trials for binary distribution

  # Initialize training data
  A <- matrix(sample(c(-1,1), size = n*1, replace = TRUE), nrow = n, ncol = 1)

  epsH <- matrix(stats::rnorm(n=n*1, mean=0, sd=1), nrow = n, ncol = 1)
  H <- epsH

  epsX <- matrix(stats::rnorm(n=n*2, mean=0, sd=1), nrow = n, ncol = 1)
  X <- matrix(nrow = n, ncol = 1)
  X <- g1*A+H+epsX

  Y <- matrix(stats::rpois(n=n, exp(X+H)), nrow = n, ncol = 1)

  # Set up test bench
  P.A <- A%*%solve(t(A)%*%A)%*%t(A) # TO WITH LM

  AGLM <- function(xi){

    # Step 1. Define the objective loss
    loss <- function(b.hat){
      return(sum(Y*(X*b.hat)-exp(X*b.hat)))
    }

    # Step 2. Define anchor penalty
    anchor_penalty <- function(b.hat){
      mu.hat <- exp(X*b.hat) # inverse of logit link

      special.case <- function(Y,mu.hat){
        ifelse(Y==0|mu.hat==0, mu.hat, Y*log(Y/mu.hat)-(Y-mu.hat))
      }

      r.D <- sign(Y-mu.hat)*sqrt(2*special.case(Y,mu.hat))
      fit.temp <- lm(r.D~A)
      return(sum((fitted(fit.temp))^2))
      #return(t(r.D)%*%P.A%*%r.D)
    }

    # Step 3. Contruct objective by 1. and 2.
    objective <- function(b.hat){
      return(1/n*(-loss(b.hat) + xi * anchor_penalty(b.hat)))
    }

    # Set start value for optimization
    ans2 <- optim(f=objective, par = runif(1), method = "L-BFGS-B")
    return(ans2$par)
  }

  xi=0
  fit.glm <- glm(formula = Y~X-1, family = poisson)
  fit.aglm <- anchorglm(formula = Y~X-1, A_formula = ~A-1, xi=xi, m=1, family=poisson, type="deviance")
  logLik(fit.glm)
  logLik(fit.aglm)

  xi = 2
  AGLM(xi)
  as.numeric(anchorglm(formula = Y~X-1, A_formula = ~A-1, xi=xi, m=1, family=poisson, type="deviance")$optim$par)

  expect_equal(AGLM(xi=2),
               as.numeric(anchorglm(formula = Y~X-1, A_formula = ~A-1, xi=xi, m=1, family=poisson, type="deviance")$optim$par), tolerance = 0.00001)

})

#####################################################################################
test_that("Testing construction of normal anchor objective and optimization", {

  set.seed(1)

  n <- 1000 # number of samples from unpertubed and pertubed distribution

  # Initialize data

  # Anchor coefficients
  g1 <- 1

  # Initialize training data
  A <- matrix(sample(c(-1,1), size = n, replace = TRUE), nrow = n, ncol = 1)

  epsH <- matrix(stats::rnorm(n=n*1, mean=0, sd=1), nrow = n, ncol = 1)
  H <- epsH

  epsX <- matrix(stats::rnorm(n=n, mean=0, sd=1), nrow = n, ncol = 1)
  X <- matrix(nrow = n, ncol = 1)
  X[,1] <- g1*A+H+epsX

  epsY <- epsY <- matrix(stats::rnorm(n=n, mean=0, sd=1), nrow = n, ncol = 1)
  Y <- matrix(X+2*H+epsY, nrow = n, ncol = 1)

  # Set up test bench
  anchor.regression <- function(X, Y, A, gamma, n){

    P.A <- A%*%solve(t(A)%*%A)%*%t(A)
    W <- diag(n)-(1-sqrt(gamma))*P.A

    Y.tilde <- W%*%Y
    X.tilde <- W%*%X

    fit <- lm(Y.tilde~X.tilde-1)
  }

  xi=0
  fit.glm <- glm(formula = Y~X-1, family = gaussian)
  fit.aglm <- anchorglm(formula = Y~X-1, A_formula = ~A-1, xi=xi, m=1, family=gaussian, type="deviance")
  logLik(fit.glm)
  logLik(fit.aglm)

  xi = 2
  gamma <- xi+1
  as.numeric(anchor.regression(X, Y, A, gamma, n)$coef)
  as.numeric(anchorglm(formula = Y~X-1, A_formula = ~A-1, xi=xi, type="deviance")$optim$par)

  # Compare results
  expect_equal(as.numeric(anchor.regression(X, Y, A, gamma, n)$coef),
               as.numeric(anchorglm(formula = Y~X-1, A_formula = ~A-1, xi=xi, type="deviance")$optim$par), tolerance = 0.01)
})

