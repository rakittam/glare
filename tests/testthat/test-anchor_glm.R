test_that("workplace for testing main anchor function", {

  set.seed(1)

  n <- 1000 # number of samples from unpertubed and pertubed distribution

  # Initialize data

  # Anchor coefficients
  g1 <- 0.5
  g2 <- -0.2
  g3 <- -2

  # initialize training data
  A <- matrix(nrow = n, ncol = 2)
  H <- matrix(nrow = n, ncol = 1)
  X <- matrix(nrow = n, ncol = 2)
  Y <- matrix(nrow = n, ncol = 1)

  m <- 5 # number of trials for binary distribution

  for (i in 1:n) {

    A[i,] <- sample(c(-1,1), size = 2, replace = TRUE) # rademacher

    epsH <- stats::rnorm(n=1, mean=0, sd=1)
    H[i] <- epsH

    epsX <- stats::rnorm(n=2, mean=0, sd=1)
    X[i,] <- g1*A[i,1]+g2*A[i,2]+H[i]+epsX

    Y[i] <- stats::rbinom(n=1, size=m, stats::plogis(3*X[i,1]+3*X[i,2]+H[i]+g3*A[i,1]))
  }

  P.A <- A%*%solve(t(A)%*%A)%*%t(A)

  AGLM <- function(xi){

    # Step 1. Define the objective loss
    loss <- function(b.hat){
      return(-sum(Y*(X%*%b.hat)-m*log(1+exp(X%*%b.hat))))
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

      #r.D <- (Y/m-p.hat)/abs(Y/m-p.hat)*sqrt(2*(Y*log(Y/(m*p.hat))+(m-Y)*log((m-Y)/(m-m*p.hat)))) # deviance residuals
      r.D <- sign(Y/m-p.hat)*sqrt(2*(special.case1(Y)+special.case2(Y))) # deviance residuals
      return(t(r.D)%*%P.A%*%r.D)
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

  expect_equal(AGLM(2), as.numeric(anchor_glm(Y=Y, X=X, A=A, xi=2, m=m, family=binomial, type="deviance")$optim$par), tolerance = 0.0001)







#
#   epsH <- rnorm(n)
#   epsX <- rnorm(n)
#   epsY <- rnorm(n)
#
#   A <- sample(c(-1,1), size = n, replace = TRUE)
#   H <- epsH
#   X <- A+H+epsX
#   Y <- X+2*H+epsY
#
#   X <- as.matrix(X)
#
#   anchor.regression <- function(X, Y, A, gamma, n){
#
#     P.A <- A%*%solve(t(A)%*%A)%*%t(A)
#     W <- diag(n)-(1-sqrt(gamma))*P.A
#
#     Y.tilde <- W%*%Y
#     X.tilde <- W%*%X
#
#     fit <- lm(Y.tilde~X.tilde-1)
#   }
#
#   expect_equal(as.numeric(anchor.regression(X, Y, A, 2, n)$coef), as.numeric(anchor_glm(Y=Y, X=X, A=A, xi=2, family=gaussian)$par), tolerance = 0.0001)
#









  # Anchor coefficients
  g1 <- 0.5

  # initialize training data
  A.train <- matrix(nrow = n, ncol = 1)
  H.train <- matrix(nrow = n, ncol = 1)
  X.train <- matrix(nrow = n, ncol = 10)
  Y.train <- matrix(nrow = n, ncol = 1)

  m <- 5 # number of trials for binary distribution

  for (i in 1:n) {

    A.train[i] <- sample(c(-1,1), size = 1, replace = TRUE)

    epsH.train <- rnorm(n=1, mean=0, sd=1)
    H.train[i] <- epsH.train
    epsX.train <- rnorm(n=10, mean=0, sd=1)
    X.train[i,] <- g1*A.train[i]+H.train[i]+epsX.train

    Y.train[i] <- rpois(n=1, exp(1*X.train[i,2]+H.train[i]))
  }

  # Objective data
  X <- X.train[,2]
  Y <- Y.train
  H <- H.train
  A <- A.train

  # Orthogonal projection on column space of anchor A
  P.A <- A%*%solve(t(A)%*%A)%*%t(A)

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
      return(t(r.D)%*%P.A%*%r.D)
    }

    # Step 3. Contruct objective by 1. and 2.
    objective <- function(b.hat){
      return(1/n*(-loss(b.hat) + xi * anchor_penalty(b.hat)))
    }

    # Set start value for optimization
    ans2 <- optim(f=objective, par = runif(1), method = "L-BFGS-B")
    return(ans2$par)
  }

  X <- as.matrix(X)
  expect_equal(AGLM(2), as.numeric(anchor_glm(Y=Y, X=X, A=A, xi=2, m=m, family=poisson, type="deviance")$optim$par), tolerance = 0.0001)

})
