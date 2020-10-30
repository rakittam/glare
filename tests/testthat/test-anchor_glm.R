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

  expect_equal(AGLM(2), anchor_glm(Y, X, A, 2, m, "binomial")$par, tolerance = 0.0001)

})
