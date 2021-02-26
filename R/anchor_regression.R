#' Anchor Regression
#'
#' Anchor regression function taken from Dominik Rothenhaeusler, Nicolai
#'  Meinshausen, Peter Buehlmann and Jonas Peters, 2021 - "Anchor regression:
#'  heterogeneous data meet causality", retrieved February 24, 2021,
#'  [https://doi.org/10.1111/rssb.12398].
#'
#' @param formula an object of class \code{"\link{formula}"} for the response
#'  and covariate variables.
#' @param A_formula an object of class \code{"\link{formula}"} for the anchor
#'  variables.
#' @param data an optional data frame, list or environment containing the
#'  variables in the model. If not found in data, the variables are taken from
#'   `environment(formula)`.
#' @param gamma numeric value for hyperparameter gamma.
#' @param ... further arguments from other methods.
#'
#' @details We recommand to use `anchor_regression` for a classical linear
#'  gaussian setup. However, there exists a relationship to `glare` using
#'  deviance residuals. The relationship and an example can be found below.
#'
#' @return returns fitted lm object of the transformed variables.
#' @export
#' @examples
#'  # Generate Data
#'  data_list <- generate_example_data(n = 1000, dim_X = 1, dim_A = 1, dim_H = 1,
#'                                     family = gaussian, A_distr = "rademacher",
#'                                     XH = 1, XA = 1, GX = 1, GH = 2)
#'  data <- data.frame(Y = data_list$Y, X = data_list$X, A = data_list$A)
#'
#'  # Set hyperparameter relation
#'  xi <- 10
#'  gamma <- 2 * xi + 1
#'
#'  # Gives approximately the same result
#'  glare(formula = Y ~ X - 1, A_formula = ~ A - 1, data = data,
#'        xi = xi, type = "deviance")$coefficients
#'
#'  anchor_regression(formula = Y ~ X - 1, A_formula = ~ A - 1, data = data,
#'                    gamma = gamma)$coefficients
anchor_regression <- function(formula, A_formula, data, gamma, ...){

  mf <- model.frame(formula, data = data)
  Y <- model.response(mf)
  X <- model.matrix(formula, data = data)
  A <- model.matrix(A_formula, data)

  n <- length(Y)

  P.A <- A%*%solve(t(A)%*%A)%*%t(A)
  W <- diag(n)-(1-sqrt(gamma))*P.A

  Y.tilde <- W%*%Y
  X.tilde <- W%*%X

  fit <- lm(Y.tilde~X.tilde-1)
}
