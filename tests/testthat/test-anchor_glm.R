test_that("workplace for testing main anchor function", {

  set.seed(1)

  n <- 1000 # number of samples from unpertubed and pertubed distribution

  # Initialize data
  library(extraDistr) # for rademacher distribution

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

    A[i,] <- extraDistr::rsign(n=2)

    epsH <- rnorm(n=1, mean=0, sd=1)
    H[i] <- epsH

    epsX <- rnorm(n=2, mean=0, sd=1)
    X[i,] <- g1*A[i,1]+g2*A[i,2]+H[i]+epsX

    Y[i] <- rbinom(n=1, size=m, plogis(3*X[i,1]+3*X[i,2]+H[i]+g3*A[i,1]))
  }

  anchor_glm(Y, X, A, 2, m, "binomial")

})
