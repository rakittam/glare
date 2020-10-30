aglm_testdata <- function(b, Y, X, m, family){

  n <- length(Y)

  test_loglikelihood <- switch(family,
                               binary = 1/n*sum(Y*(X%*%b)-m*log(1+exp(X%*%b))),
                               poisson = 1/n*sum(Y*(X%*%b)-exp(X%*%b))
                               )

  return(test_loglikelihood)
}

# besser zweimal die binary_likelihood etc. functions benutzen mit passenden argumenten
