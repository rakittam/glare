test_that("Testing binomial anchor objective and optimization", {

  set.seed(1)

  # Generate example data -----------------------------------------------------

  # H = epsH
  HX = NULL
  HA = NULL

  # X1 = H + 0.5 * A1 - 0.2 * A2 + epsX
  # X2 = H + 0.5 * A1 - 0.4 * A2 + epsX
  XH = rbind(1, 1)
  XA = rbind(c(0.5, -0.2), c(0.5, -0.4))

  # G = 3 * X1 + 3 * X2 + H - 2 * A1
  GX = c(3, 3)
  GH = 1
  GA = c(-2, 0)

  data_list <- generate_example_data(n = 1000, dim_X = 2, dim_A = 2, dim_H = 1,
                                     family = binomial, m_max = 5,
                                     A_distr = "rademacher",
                                     XH = XH, XA = XA,
                                     HX = HX,  HA = HA,
                                     GX = GX, GH = GH, GA = GA)

  data <- data.frame(Y = data_list$Y,
                     X = data_list$X,
                     A = data_list$A,
                     m = data_list$m)

  # Test bench ----------------------------------------------------------------

  # Set up data
  Y <-  data_list$Y
  X <-  data_list$X
  A <-  data_list$A
  m <-  data_list$m
  n <- length(Y)

  # Set up test bench function
  AGLM <- function(xi) {

    # Step 1. Define the objective loss
    loss <- function(b.hat) {
      -sum(log(choose(m, Y)) +
             Y * (X %*% b.hat) - m * log(1 + exp(X %*% b.hat)))
    }

    # Step 2. Define anchor penalty
    anchor_penalty <- function(b.hat) {
      # inverse of logit link
      p.hat <- exp(X %*% b.hat) / (1 + exp(X %*% b.hat))

      special.case1 <- function(Y) {
        ifelse(Y == 0, 0, Y * log(Y / (m * p.hat)))
      }
      special.case2 <- function(Y) {
        ifelse(Y == m, 0, (m - Y) * log((m - Y) / (m - m * p.hat)))
      }

      # deviance residuals
      r.D <- sign(Y / m - p.hat) *
        sqrt(2 * (special.case1(Y) + special.case2(Y)))

      fit_temp <- lm(r.D ~ A)
      sum((fitted(fit_temp))^2)
    }

    # Step 3. Construct objective by 1. and 2.
    objective <- function(b.hat) {
      1 / n * (loss(b.hat) + xi * anchor_penalty(b.hat))
    }

    # Step 4. Optimization
    ans2 <- optim(fn = objective, par = runif(2), method = "L-BFGS-B")
    ans2$par
  }

  # Compare methods -----------------------------------------------------------

  # log-likelihood test
  xi <-  0
  YY <- cbind(Y, m - Y)
  fit_glm <- glm(formula = YY ~ X - 1, family = binomial)
  fit_glare <- glare(formula = YY ~ X.1 + X.2 - 1,
                     A_formula = ~ A.1 + A.2 - 1,
                     data = data,
                     xi = xi,
                     family = binomial,
                     type = "deviance")

  logLik_glm <- as.numeric(logLik(fit_glm))
  logLik_glare <- logLik(fit_glare)
  logLik_glm
  logLik_glare

  expect_equal(logLik_glm, logLik_glare, tolerance = 0.001)

  # Parameter estimation test
  xi <-  2
  b_AGLM <- AGLM(xi = xi)
  b_glare <- glare(formula = Y ~ X.1 + X.2 - 1,
                   A_formula = ~ A.1 + A.2 - 1,
                   data = data,
                   xi = xi,
                   family = binomial,
                   type = "deviance")$optim$par
  b_AGLM
  b_glare

  expect_equal(b_AGLM, b_glare, tolerance = 0.001)

  # (success, fails) input test
  YY <- cbind(Y, m-Y)
  b_glare_YY <- glare(formula = YY ~ X.1 + X.2 - 1,
                      A_formula = ~ A.1 + A.2 - 1,
                      data = data,
                      xi = xi,
                      family = binomial,
                      type = "deviance")$optim$par
  b_glare_YY

  expect_equal(b_glare, b_glare_YY, tolerance = 0.001)
})

###############################################################################
test_that("Testing poisson anchor objective and optimization", {

  set.seed(1)

  # Generate example data -----------------------------------------------------

  # H = epsH
  HX = NULL
  HA = NULL

  # X = H + 0.5 * A + epsX
  XH = 1
  XA = 0.5

  # G = X + H
  GX = 1
  GH = 1
  GA = NULL

  data_list <- generate_example_data(n = 1000, dim_X = 1, dim_A = 1, dim_H = 1,
                                     family = poisson,
                                     A_distr = "rademacher",
                                     XH = XH, XA = XA,
                                     HX = HX,  HA = HA,
                                     GX = GX, GH = GH, GA = GA)

  data <- data.frame(Y = data_list$Y,
                     X = data_list$X,
                     A = data_list$A)

  # Test bench ----------------------------------------------------------------

  # Set up data
  Y <-  data_list$Y
  X <-  data_list$X
  A <-  data_list$A
  n <- length(Y)

  # Set up test bench function
  AGLM <- function(xi) {

    # Step 1. Define the objective loss
    loss <- function(b.hat){
      sum(Y * (X * b.hat) - exp(X * b.hat))
    }

    # Step 2. Define anchor penalty
    anchor_penalty <- function(b.hat) {
      mu.hat <- exp(X * b.hat) # inverse of logit link

      special.case <- function(Y,mu.hat) {
        ifelse(Y == 0 | mu.hat == 0, mu.hat, Y * log(Y / mu.hat) - (Y - mu.hat))
      }

      r.D <- sign(Y - mu.hat) * sqrt(2 * special.case(Y, mu.hat))
      fit_temp <- lm(r.D ~ A)
      sum((fitted(fit_temp))^2)
    }

    # Step 3. Construct objective by 1. and 2.
    objective <- function(b.hat) {
      1 / n * (-loss(b.hat) + xi * anchor_penalty(b.hat))
    }

    # Step 4. Optimization
    ans2 <- optim(fn = objective, par = runif(1), method = "L-BFGS-B")
    ans2$par
  }

  # Compare methods -----------------------------------------------------------

  # log-likelihood test
  xi <-  0
  fit_glm <- glm(formula = Y ~ X - 1, family = poisson)
  fit_glare <- glare(formula = Y ~ X - 1,
                     A_formula = ~ A - 1,
                     data = data,
                     xi = xi,
                     family = poisson,
                     type = "deviance")

  logLik_glm <- as.numeric(logLik(fit_glm))
  logLik_glare <- logLik(fit_glare)
  logLik_glm
  logLik_glare

  expect_equal(logLik_glm, logLik_glare, tolerance = 0.001)

  # Parameter estimation test
  xi <-  2.5
  b_AGLM <- AGLM(xi = xi)
  b_glare <- glare(formula = Y ~ X - 1,
                   A_formula = ~ A - 1,
                   data = data,
                   xi = xi,
                   family = poisson,
                   type = "deviance")$optim$par
  b_AGLM
  b_glare

  expect_equal(b_AGLM, b_glare, tolerance = 0.001)
})

###############################################################################
test_that("Testing construction of normal anchor objective and optimization", {

  set.seed(1)

  # Generate example data -----------------------------------------------------

  # H = epsH
  HX = NULL
  HA = NULL

  # X = H + A + epsX
  XH = 1
  XA = 1
  # G = X + 2 * H (Y = G + epsY)
  GX = 1
  GH = 2
  GA = NULL

  data_list <- generate_example_data(n = 1000, dim_X = 1, dim_A = 1, dim_H = 1,
                                     family = gaussian,
                                     A_distr = "rademacher",
                                     XH = XH, XA = XA,
                                     HX = HX,  HA = HA,
                                     GX = GX, GH = GH, GA = GA)

  data <- data.frame(Y = data_list$Y,
                     X = data_list$X,
                     A = data_list$A)

  # Test bench ----------------------------------------------------------------

  # Set up data
  Y <-  data_list$Y
  X <-  data_list$X
  A <-  data_list$A

  n <- length(Y)

  # Set up test bench function
  anchor.regression <- function(X, Y, A, gamma, n) {

    P.A <- A %*% solve(t(A) %*% A) %*% t(A)
    W <- diag(n) - (1 - sqrt(gamma)) * P.A

    Y.tilde <- W %*% Y
    X.tilde <- W %*% X

    fit <- lm(Y.tilde ~ X.tilde - 1)
    as.numeric(fit$coefficients)
  }

  # Compare methods -----------------------------------------------------------

  # log-likelihood test
  xi <-  0
  fit_glm <- glm(formula = Y ~ X - 1, family = gaussian)
  fit_glare <- glare(formula = Y ~ X - 1,
                     A_formula = ~ A - 1,
                     data = data,
                     xi = xi,
                     family = gaussian,
                     type = "deviance")

  logLik_glm <- as.numeric(logLik(fit_glm))
  logLik_glare <- logLik(fit_glare)

  logLik_glm
  logLik_glare

  expect_equal(logLik_glm, logLik_glare, tolerance = 0.001)

  # Parameter estimation test
  xi <- 0.5
  b_AGLM <- anchor.regression(X, Y, A, gamma = xi + 1, n)
  b_glare <- glare(formula = Y ~ X - 1,
                   A_formula = ~ A - 1,
                   data = data,
                   xi = xi,
                   type = "deviance")$optim$par

  b_AGLM
  b_glare

  expect_equal(b_AGLM, b_glare, tolerance = 0.001)

  # IV test
  xi <- 1000000
  b_AGLIV <- anchor.regression(X, Y, A, gamma = xi + 1, n)
  b_glare_IV <- glare(formula = Y ~ X - 1,
                      A_formula = ~ A - 1,
                      data = data,
                      xi = xi,
                      type = "deviance")$optim$par

  b_AGLIV
  b_glare_IV

  expect_equal(b_AGLIV, b_glare_IV, tolerance = 0.001)
})
