Generalized Linear Anchor Regression
================
Maic Rakitta
2020-12-05

<!-- README.md is generated from README.Rmd. Please edit that file -->

This vignette offers an introduction with examples to the first
Generalized Linear Anchor REgression (GLARE) package version. We
extended the work of Dominik Rothenhaeusler, Nicolai Meinshausen, Peter
Buehlmann and Jonas Peters, who developed a linear anchor regression
method. GLARE extends their work from a linear to a generalized linear
framework. In this vignette we give only a thin theoretical introduction
to GLARE and focus more on how the reader can apply the corresponding R
package. Since the development is not yet finished, we are grateful for
your input. Furthermore, we cannot yet provide GLARE literature. For
anchor regression we recommend the corresponding work of Rothenhaeusler
et al.

This vignette makes no claim to completeness, but rather provides a
first insight into our method.

## Short Introduction

Similarly to the theory of Rothenhaeusler et al., GLARE estimates its
parameter by minimizing an objective function. However, instead of
interpolating between the solutions of ordinary least squares and
two-stage least squares, we adapt GLARE to the GLM framework. The main
differences are that we use log-likelihood instead of OLS and that
suitable residuals must be used. The minimization results in
\[\beta^\xi = \arg\min_\beta -\frac{1}{n}l(\beta; Y, X) +
\frac{1}{n}\xi ||Ar||_2^2,\] where \(X \in \mathbb{R}^{n\times p}\) are
the covariates, \(Y \in \mathbb{R}^{n}\) is the response,
\(A \in \mathbb{R}^{n\times q}\) are the anchor variables and
\(r \in \mathbb{R}^{n}\) are appropriate residuals (e.g. deviance
residuals, pearson residuals). The hyperparameter
\(\xi \in \{-1, \infty\}\) interpolates between maximizing the
log-likelihood, \(l\), and partialling out the effect of the anchor onto
the residuals. Note that there can also be hidden variables,
\(H \in \mathbb{R}^{n\times h}\), influencing both the response and the
covariates.

## Getting started

The package is still in development. For the most up-to-date version,
feel free to send a request directly to the author
(<rakittam@student.ethz.ch>) to receive the corresponding source
package.

If you need help to install the package from source, please have a look
at <http://www.ryantmoore.org/files/ht/htrtargz.pdf> .

For this vignette, up to GLARE, the following libraries are necessary:

``` r
library(glare) 
library(stats) # used to generate data sets
library(ggplot2) # used for plots
```

## Example: Binomial Framework

The given example is constructed in a binomial framework. First we
construct an artificial data set. Next, we show how to apply GLARE to
the generated data. Finally, we investigate the model on a perturbed
data set and plot some result.

### Data Construction

In this section we simulate the unperturbed data for the binomial
example. The chosen model for this example is \[
\begin{aligned}
H &\sim \mathcal{N}(0, 1)\\
A &\sim \text{Rademacher}_2\\
X &= M_1 A + H + \epsilon_X, \text{ where } \epsilon_X \sim \mathcal{N}_2(0, 1)\\
Y&\sim \text{binomial}(p), \text{ where } g(p) = \text{logit}(p) = b^TX + H + M_2 A\text{ .} \end{aligned}
\]

If you want to reproduce our results, please set the corresponding seed,
`set.seed(1992)`, before you generate the data. Please have a look at
the help file , `?generate_example_data`, for detailed description of
the needed input.

``` r
set.seed(1992)

# Generate unperturbed data ---------------------------------------------------

  # Number of observations
  n <- 1000 

  # H = epsH
  HX = NULL
  HA = NULL

  # X1 = H + 0.5 * A1 - 0.2 * A2 + epsX
  # X2 = H + 0.5 * A1 - 0.4 * A2 + epsX
  XH = rbind(1, 1)
  XA = rbind(c(0.5, -0.2), c(0.5, -0.4))

  # G = 0.4 * X1 + 0.6 * X2 + 0.4 * H - 0.8 * A1
  GX = c(0.4, 0.6)
  GH = 0.4
  GA = c(-0.8, 0)

  data_list <- generate_example_data(n = n, dim_X = 2, dim_A = 2, dim_H = 1,
                                     family = binomial,
                                     A_distr = "rademacher",
                                     XH = XH, XA = XA,
                                     HX = HX,  HA = HA,
                                     GX = GX, GH = GH, GA = GA)

  data <- data.frame(Y = data_list$Y,
                     X = data_list$X,
                     A = data_list$A)
```

### How to apply `glare`

Now that we have an artificial training data, we can apply GLARE.
Response and covariate variables are fed into the method using the
`formula` class. Also the anchors are provided with a corresponding
`formula`. Note that we are not using an intercept, and hence added
`- 1` to the formula expression. We encourage the reader to have a look
at the help file, `?glare`, for further information.

For a first intuition we set the hyperparameter, \(\xi\), to two and use
deviance residuals for the anchor penalty.

``` r
glare_fit <- glare(formula = Y ~ X.1 + X.2 - 1,
                  A_formula = ~ A.1 + A.2 - 1,
                  data = data,
                  xi = 2,
                  family = "binomial",
                  type = "deviance")
```

By specifying an anchor variable, a symbolic description of the linear
predictor and a description of the error distribution, we have fitted a
generalized linear anchor model. The function returns an object of class
`"glare"`.

It shall be noted that `glare` can handle different form of input for
binomial data, i.e.

``` r
# Response as a matrix with number of success and losses for each observation
SL <- cbind(data_list$Y, 1 - data_list$Y)
data_SL <- list(SL = SL, X = data_list$X, A = data_list$A)

glare_fit_SL <- glare(formula = SL ~ X - 1,
                     A_formula = ~ A - 1,
                     data = data_SL,
                     xi = 2,
                     family = "binomial",
                     type = "deviance")
```

### Interpretation on a perturbed data set

In this section we want to investigate our model trained on the initial
data using different values for \(\xi\), on a perturbed data set.

#### Perturbed data set

First we simulate a perturbed data set by manipulating only the anchor
coefficients, while letting the rest the same as in the initial data. To
see the changes compared to the original undisturbed data set, we write
down again the formulas from above that change. \[
\begin{aligned}
X_1 &= H + 0.5A_1 - 0.2A_2 + \epsilon_{X_1}\\
X_1^\text{pert} &= H - 1.5A_1 - 0.5A_2 + \epsilon_{X_1^\text{pert}}\\
X_2 &= H + 0.5A_1 - 0.4A_2 + \epsilon_{X_2}\\
X_2^\text{pert} &= H - 1.5A_1 - 0.4A_2 + \epsilon_{X_2^\text{pert}}
\end{aligned}
\begin{aligned}
G &= 3X_1 + 3X_2 + H - 2A_1\\
G^\text{pert} &= 3X_1^\text{pert} + 3X_2^\text{pert} + H + 2A_1
\end{aligned}
\]

``` r
# Generate perturbed data ---------------------------------------------------

# X1_pert = H - 1.5 * A1 - 0.5 * A2 + epsX
# X2_pert = H - 1.5 * A1 - 0.4 * A2 + epsX
XA_pert = rbind(c(-1.5, -0.5), c(-1.5, -0.4))

# G_pert = 0.4 * X1 + 0.6 * X2 + 0.4 * H + 0.8 * A1
GA_pert = c(0.8, 0)

data_list_pert <- generate_example_data(n = n, dim_X = 2, dim_A = 2, dim_H = 1,
                                        family = binomial,
                                        A_distr = "rademacher",
                                        XH = XH, XA = XA_pert,
                                        HX = HX,  HA = HA,
                                        GX = GX, GH = GH, GA = GA_pert)

data_pert <- data.frame(Y = data_list_pert$Y,
                        X = data_list_pert$X,
                        A = data_list_pert$A)
```

#### Iterating over the hyperparameter

Next, we iterate over \(\xi\) and use `glare` for each iteration. We
store the estimated coefficients and the log-likelihood on the perturbed
data set.

``` r
xi_vec <- seq(-1, 5, by = 0.1)
b_matrix <- matrix(nrow = length(xi_vec), ncol = 2)
loglikelihood_pert <- numeric(length(xi_vec))

for (i in 1:length(xi_vec)) {

  xi <- xi_vec[i]
  fit_temp <- glare(formula = Y ~ X.1 + X.2 - 1,
                    A_formula = ~ A.1 + A.2 - 1,
                    data = data,
                    xi = xi,
                    family = "binomial",
                    type = "deviance")

  b_matrix[i, ] <- coef(fit_temp)
  loglikelihood_pert[i] <- logLik(fit_temp, newdata = data_pert)
}
```

Moreover, we calculate additional estimates for comparison: maximum
likelihood estimation, an analogy to partialling out the effect of the
anchors A and letting \(\xi\) go to infinity. In the following the
corresponding parameter estimates and the corresponding log-likelihood
on the perturbed data is given. The reader should be aware that like
anchor regression, GLARE does not aim for the true underlying causal
function, but rather tries to maximize the log-likelihood.

``` r
# Optimal xi ------------------------------------------------------------------
xi_opt <- xi_vec[which.max(loglikelihood_pert)]
xi_opt
#> [1] 1.2
b_opt <- b_matrix[which.max(loglikelihood_pert), ]
b_opt
#> [1] 0.2825194 0.4221096
loglikelihood_pert_opt <- max(loglikelihood_pert)
loglikelihood_pert_opt
#> [1] -523.1204

# MLE -------------------------------------------------------------------------
xi_MLE <- 0
b_MLE <- b_matrix[which(xi_vec == xi_MLE), ]
b_MLE
#> [1] 0.4325611 0.5144322
# which by construction is the same as
glm(SL ~ data_list$X - 1, family = "binomial")$coef
#> data_list$X1 data_list$X2 
#>    0.4325611    0.5144322
loglikelihood_pert_MLE <- loglikelihood_pert[which(xi_vec == xi_MLE)]
loglikelihood_pert_MLE 
#> [1] -535.1648

# PA-like ---------------------------------------------------------------------
xi_PA <- -1
b_PA <- b_matrix[which(xi_vec == xi_PA), ]
b_PA
#> [1] 0.6350117 0.6979080
loglikelihood_pert_PA <- loglikelihood_pert[which(xi_vec == xi_PA)]
loglikelihood_pert_PA
#> [1] -588.5285

# big xi ----------------------------------------------------------------------
xi_big <- 1000
fit_big <- glare(formula = Y ~ X.1 + X.2 - 1,
                 A_formula = ~ A.1 + A.2 - 1,
                 data = data,
                 xi = xi_big,
                 family = "binomial",
                 type = "deviance")
b_big <- coef(fit_big)
b_big
#>       X.1       X.2 
#> -1.594340  1.586996
loglikelihood_pert_big <- logLik(fit_big, newdata = data_pert)
loglikelihood_pert_big
#> [1] -1119.938
```

#### Plot of the results

Finally, we can plot the log-likelihood on the perturbed data against
the hyperparameter values. As we can see in the figure, in this setup
GLARE outperforms the other methods with the right hyperparameter,
\(\xi\), in terms of maximal log-likelihood on the perturbed data.

``` r
ggplot_data <- data.frame(loglikelihood_pert = loglikelihood_pert,
                          xi_vec = xi_vec)
ggplot(data = ggplot_data, aes(y = loglikelihood_pert, x = xi_vec)) +

  geom_line() +
  
  labs(title = "log-likelihood of the perturped data set",
       x = "Hyperparameter xi", y = "log-likelihood") +

  annotate("point", colour = "red",
           x = xi_opt, y = loglikelihood_pert_opt) +
  annotate("text", label = "optimal xi",
           x = xi_opt, y = loglikelihood_pert_opt + 1.5) +
  annotate("point", colour = "green4",
           x = xi_MLE, y = loglikelihood_pert_MLE) +
  annotate("text", label = "MLE",
           x = xi_MLE + 0.3, y = loglikelihood_pert_MLE) +
  annotate("point", colour = "blue1",
           x = xi_PA, y = loglikelihood_pert_PA) +
  annotate("text", label = "PA-like",
           x = xi_PA + 0.4, y = loglikelihood_pert_PA) 
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

### Methods

We provide several generic accessor functions for the `glare` class,
i.e. `logLik`, `coef`, `predict` and `residuals`. For a list of all
currently provided methods have a look at `methods(class = "glare")`.

For example you can generate predictions given by the learned model for
either the initial or for new data on the link or the response scale.

``` r
predict(glare_fit, type = "link") # predictions on link scale
predict(glare_fit, type = "response") # predictions on response scale
predict(glare_fit, newdata = data_pert, type = "link") # perturbed data
```

Or one can print the summary of a fitted `glare` object.

``` r
summary(glare_fit)
#> 
#> Call:  glare(formula = Y ~ X.1 + X.2 - 1, A_formula = ~A.1 + A.2 - 1, 
#>     data = data, xi = 2, family = "binomial", type = "deviance")
#> 
#> Deviance Residuals:
#>      Min.  1st Qu.  Median    Mean  3rd Qu.   Max.
#>   -1.8938  -0.9899  0.5964  0.0487   1.0051  1.951
#> 
#> Coefficients:
#>      Estimate  Std. Error  t value      Pr(>|t|)
#> X.1    0.2105      0.0514   4.0957  4.209421e-05
#> X.2    0.3966      0.0507   7.8204  5.263683e-15
```
