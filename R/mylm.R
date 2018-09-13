
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
  est$data <- X

  #Fit model (using LSQ)
  est$coefficients = solve((t(X) %*% X), t(X) %*% y)

  #Construct fitted values
  est$fitted_values = X%*%est$coefficients
  est$residuals = y - est$fitted_values

  #Estimate error
  est$SSE = t(est$residuals) %*% est$residuals
  est$estimated_variance = est$SSE / (dim(X)[1] - dim(X)[2])

  #store covariance matrix of parameter estimates
  est$cov_coeff = solve(t(X)%*%X) * as.numeric(est$estimated_variance)


  #Compute test statistics
  est$std_coeff = sqrt(diag(est$cov_coeff))
  est$coeff_z = est$coefficients/est$std_coeff
  est$p_values = pmax(2*(1-pnorm(abs(est$coeff_z))),rep(2e-16,dim(X)[2]))

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
}

plot.mylm <- function(object, ...){
  # Code here is used when plot(object) is used on objects of class "mylm"
  cat("Funny pictures")
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
