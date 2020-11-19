#' Generate example data for glare
#'
#' @param n number of observations.
#' @param dim_X number of covariates.
#' @param dim_A number of anchor variables.
#' @param dim_H number of hidden confounders.
#' @param family a description of the error distribution and link function to
#'  be used in the model. This can be a character string naming a family
#'  function, a family function or the result of a call to a family function.
#'  See \code{\link{family}} for details. The families supported up to now are
#'  `gaussian`, `binomial` and `poisson`.
#' @param m_max maximal number of indepdent bernoulli trials, if binomial family
#'  is used.
#' @param A_distr distribution from which the Anchor variable is generated.
#' @param XH coefficients of how X depends on H as a dim_X times dim_H matrix.
#' @param XA coefficients of how X depends on A as a dim_X times dim_A matrix.
#' @param HX coefficients of how H depends on X as a dim_H times dim_X matrix.
#' @param HA coefficients of how H depends on A as a dim_H times dim_A matrix.
#' @param GX coefficients of how G depends on X as a dim_X vector.
#' @param GH coefficients of how G depends on H as a dim_H vector.
#' @param GA coefficients of how G depends on A as a dim_A vector.
#' @param ... other
#'
#' @details For this function only one direction of dependence between X and H
#'  are allowed.\cr\cr
#'  The response Y is generated with the corresponding family, linkfunction and
#'  the variable G, e.g.\cr `Y = rpois(n = n, lambda = linkinv(G))`.\cr\cr
#'  As input the desired dimension of the variables are used for a weak
#'  check for correct input of the variable coefficients.\cr\cr
#'  For a script to generate data manipulate the given example below.
#'
#' @return returns a data list of matrices of all the constructed variables.
#'
#' @examples
#'  # H = epsH
#'  HX = NULL
#'  HA = NULL
#'
#'  # X1 = H + 0.5 * A1 - 0.2 * A2 + epsX
#'  # X2 = H + 0.5 * A1 - 0.4 * A2 + epsX
#'  XH = rbind(1, 1)
#'  XA = rbind(c(0.5, -0.2), c(0.5, -0.4))
#'
#'  # G = 3 * X1 + 3 * X2 + H - 2 * A1
#'  GX = c(3, 3)
#'  GH = 1
#'  GA = c(-2, 0)
#'
#'  generate_example_data(n = 1000, dim_X = 2, dim_A = 2, dim_H = 1,
#'                        family = binomial, m_max = 5,
#'                        A_distr = "rademacher",
#'                        XH = XH, XA = XA,
#'                        HX = HX,  HA = HA,
#'                        GX = GX, GH = GH, GA = GA)
#' @export
#' @importFrom stats rnorm rbinom rpois family
generate_example_data <- function(n, dim_X = 1, dim_A = 1, dim_H = 1,
                                  family = gaussian, m_max = 1,
                                  A_distr = c("normal", "rademacher"),
                                  XH = NULL, XA = NULL,
                                  HX = NULL,  HA = NULL,
                                  GX = NULL, GH = NULL, GA = NULL, ...) {

  # Allocate glm family as in glm source code
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if (family$family == "binomial") {
    m <- sample(1:m_max, size = n, replace = TRUE)
  } else {
    m <- 1
  }

  # Assign needed functions depending on given glm family
  linkinv <- family$linkinv

  # Generate data
  A_distr <- match.arg(A_distr)
  A <- switch(A_distr,
              "normal" = matrix(rnorm(n * dim_A), nrow = n, ncol = dim_A),
              "rademacher" = matrix(sample(c(-1, 1), size = n * dim_A, replace = TRUE), nrow = n, ncol = dim_A))

  epsX <- matrix(rnorm(n * dim_X), nrow = n, ncol = dim_X)
  epsH <- matrix(rnorm(n * dim_H), nrow = n, ncol = dim_H)

  X <- matrix(NA, nrow = n, ncol = dim_X)
  H <- matrix(NA, nrow = n, ncol = dim_H)
  G <- matrix(NA, nrow = n, ncol = 1)
  Y <- matrix(NA, nrow = n, ncol = 1)

  if (is.null(XA)) XA <- matrix(0, nrow = dim_X, ncol = dim_A)
  if (is.null(HA)) HA <- matrix(0, nrow = dim_H, ncol = dim_A)
  if (is.null(GX)) GX <- matrix(0, nrow = 1, ncol = dim_X)
  if (is.null(GH)) GH <- matrix(0, nrow = 1, ncol = dim_H)
  if (is.null(GA)) GA <- matrix(0, nrow = 1, ncol = dim_A)

  if (!is.null(XH) & !is.null(HX)) {
    stop("Only one direction between X and H allowed!")
  } else if (!is.null(XH)) {
    for (i in 1:dim_H) {
      H[, i] <- as.matrix(HA)[i, ] %*% t(A) + epsH[, i]
    }
    for (i in 1:dim_X) {
      X[, i] <- as.matrix(XH)[i, ] %*% t(H) + as.matrix(XA)[i, ] %*% t(A) + epsX[, i]
    }
  } else if (!is.null(HX)) {
    for (i in 1:dim_X) {
      X[, i] <- as.matrix(XA)[i, ] %*% t(A) + epsX[, i]
    }
    for (i in 1:dim_H) {
      H[, i] <- as.matrix(HX)[i, ] %*% X + as.matrix(HA)[i, ] %*% A + epsH[, i]
    }
  }

  G <- GX %*% t(X) + GH %*% t(H) + GA %*% t(A)

  mu <- linkinv(G)

  Y <- switch(family$family,
                    "binomial" = rbinom(n = n, size = m, prob = mu),
                    "poisson" = rpois(n = n, lambda = mu),
                    "gaussian" = matrix(mu + (rnorm(n = n)), nrow = n))

  list(Y = Y, X = X, A = A, H = H, G = G, family = family, m = m)
}
