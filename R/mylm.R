
# Select Build, Build and reload to build and lode into the R-session.

mylm <- function(formula, data = list(), contrasts = NULL, ...){

  # formula = rent ~ area + etc.
  # data

  # Extract model matrix & responses
  mf <- model.frame(formula = formula, data = data)
  X  <- model.matrix(attr(mf, "terms"), data = mf, contrasts.arg = contrasts)
  y  <- model.response(mf)
  terms <- attr(mf, "terms")

  # Add code here to calculate coefficients, residuals, fitted values, etc...
  # and store the results in the list est
  est <- list(terms = terms, model = mf)

  # Store call and formula used
  est$call <- match.call()
  est$formula <- formula

  # Store information used for model fit
  est$X <- X
  est$y = y

  #Fit model (using LSQ)
  est$coefficients = solve((t(X) %*% X), t(X) %*% y)

  #Construct fitted values and residuals
  est$fitted_values = X%*%est$coefficients
  est$residuals = y - est$fitted_values

  #Construct quantiles
  est$quantiles = quantile(as.numeric(est$residuals))

  #Compute sums of squares
  est$SST = t(est$y - mean(est$y)) %*% (est$y - mean(est$y))
  est$SSR = t(est$fitted_values - mean(est$y)) %*% (est$fitted_values - mean(est$y))
  est$SSE = t(est$residuals) %*% est$residuals

  #Estimate error
  est$estimated_variance = est$SSE / (dim(X)[1] - dim(X)[2])

  #store covariance matrix of parameter estimates
  est$cov_coeff = solve(t(X)%*%X) * as.numeric(est$estimated_variance)

  #Compute T-test statistics
  est$std_coeff = sqrt(diag(est$cov_coeff))
  est$coeff_z = est$coefficients/est$std_coeff
  est$p_values = pmax(2*(1-pnorm(abs(est$coeff_z))),rep(2e-16,dim(X)[2]))
  fsign = function(value){
    if (value < 0.001) {
      return ('***')
    }
    else if (value < 0.01) {
      return ('**')
    }
    else if (value < 0.05) {
      return ('*')
    }
    else if (value < 0.1) {
      return ('.')
    }
    else {
      return (' ')
    }
  }
  est$significance_lvls = sapply(est$p_values, fsign)

  #Compute X^2 test statistics
  est$F_statistic = (est$SSR / (dim(X)[2] - 1)) / (est$estimated_variance)
  est$F_p_val = max((1-pchisq(est$F_statistic*(dim(X)[2] - 1), dim(X)[2] - 1)), 2E-16)

  # Set class name. This is very important!
  class(est) <- 'mylm'

  # Return the object with all results
  return(est)
}

print.mylm <- function(object, ...){
  # Code here is used when print(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  cat('Call: \n')
  print(object$call)
  cat('\n')
  cat('Coefficients: \n')
  t(object$coefficients)
}

summary.mylm <- function(object, ...){
  # Code here is used when summary(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  cat('Summary of object\n')
  cat('Call: \n')
  print(object$call)
  cat('\nResiduals: \n')
  print(object$quantiles)
  cat('\nCoefficients: \n')

  outmatrix = data.frame(data= c(object$coefficients, object$std_coeff, object$coeff_z, object$p_values, object$significance_lvls),
                         nrow = dim(object$X)[2], ncol=5)
  outmatrix = data.frame(format(object$coefficients, digits=4, nsmall=4),
                         format(object$std_coeff, digits=4, nsmall=4),
                         format(object$coeff_z, digits=4, nsmall=2),
                         format(object$p_values, digits=4, nsmall=2),
                         object$significance_lvls)

  #rownames(outmatrix) = rownames(A$coefficients)

  colnames(outmatrix) = c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)', ' ')
  print(outmatrix, quote = F)
  cat('Signif. codes:  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1')

  cat('\n\n')
  cat('Residual Standard Error: ')
  cat(format(sqrt(object$estimated_variance), digits=1, nsmall=1))
  cat(' on '); cat(dim(object$X)[1] - dim(object$X)[2]); cat(' degrees of freedom')

  cat('\n')
  cat('Multiple R-squared: ')
  cat(format(object$SSR / object$SST, digits = 4, nsmall= 4))
  cat('    Adjusted R-squared: ')
  cat(format((1 - (object$SSE / (dim(object$X)[1] - dim(object$X)[2]))/
                    (object$SST / (dim(object$X)[1] - 1))), digits = 5))

  cat('\n')
  cat('Chisq-statisic: ')
  cat(format(object$F_statistic, digits=1, nsmall = 0))
  cat(' on '); cat(dim(object$X)[2] - 1); cat(' degrees of freedom')
  cat(', p-value: '); cat(object$F_p_val)

  cat('\n')


}

plot.mylm <- function(object, ...){
  # Code here is used when plot(object) is used on objects of class "mylm"
  qplot(object$fitted_values, object$residuals)
}



# This part is optional! You do not have to implement anova
anova.mylm <- function(object, ...){
  # Code here is used when anova(object) is used on objects of class "mylm"

  # Components to test
  comp <- attr(object$terms, "term.labels")

  # Name of response
  response <- deparse(object$terms[[2]])

  # Fit the sequence of models
  txtFormula <- paste(response, "~", sep = "")
  model <- list()
  for(numComp in 1:length(comp)){
    if(numComp == 1){
      txtFormula <- paste(txtFormula, comp[numComp])
    }
    else{
      txtFormula <- paste(txtFormula, comp[numComp], sep = "+")
    }
    formula <- formula(txtFormula)
    model[[numComp]] <- lm(formula = formula, data = object$model)
  }

  # Print Analysis of Variance Table
  cat('Analysis of Variance Table\n')
  cat(c('Response: ', response, '\n'), sep = '')
  cat('          Df  Sum sq X2 value Pr(>X2)\n')
  for(numComp in 1:length(comp)){
    # Add code to print the line for each model tested
  }

  return(model)

}
