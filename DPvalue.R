#' Type I Error Rate and Distribution of P Values
#' @param distribution contains 8 types of distributions which are Normal Distribution, Chi-square Distribution, Exponential Distribution, LaPlace Distribution, Rayleigh Distribution, Weibull Distribution, Cauchy Distribution, and Log Normal Distribution. The particular distribution for which the data is generated from can be specified using this parameter.
#' @param simulation is the number of simulations.
#' @param n is the number of the samples.
#' @param k is an offset to ensure that all fractional ranks are strictly between zero and one.
#' @param siglevel is the significance level for the two sample test.
#' @details Use a given matrix A = [0, 0, 0 ...... 1, 1, 1] to fit the data generated from a target distribution. Then transform the data using the inverse normal transformation and fit it with the matrix A again. Calculate the type I error rates for the two tests and show the distributions of the p values to evaluate the significance of the transformation.
#' @details Need to install INT_Function package first.
#' @return a table of type I error rates of the two two sample tests conducted using original and transformed data respectively and two histograms that display the distributions of p values in the two tests.


dpvalue <- function(simulation, n, k, siglevel, distribution) {
  X1 = c(rep(1, n/2), rep(0, n/2))
  X2 = c(rep(0, n/2), rep(1, n/2))
  A = matrix(c(X1, X2), nrow = n, ncol=2) # construct the matrix A = [0, 0, 0 ...... 1, 1, 1]
  p_original = vector(length = simulation) # store the p values for original data
  error_original = NA # mark the position where type I error occurs in original data
  p_transformed = vector(length = simulation) # store the p values for transformed data
  error_transformed = NA # mark the position where type I error occurs in transformed data

  if (distribution == "norm") {
    # method of processing for all other distributions are constructed in the same way as this one
    mean = as.integer(readline(prompt = "mean : ")) # ask for the mean of the distribution
    sd = as.integer(readline(prompt = "sd : ")) # ask for the standard deviation of the distribution
    plot(density(rnorm(n, mean, sd)), main = "Density Curve") # plot the density curve for reference
    for(i in 1 : simulation){ # number of simulations
      original = rnorm(n, mean, sd) # generate data from the desired distribution
      fit_original = lm(formula = original ~ A) # perform the two sample t test using the samples generated and the given matrix A
      p_original[i] = summary(fit_original)$coefficients[2,4] # collect the p values
      if (p_original[i] < siglevel) {
        error_original[i] = 1 # mark the positions of errors in order to calculate the type I error rate
      }
      transformed = INT_function(original, n, k) # apply the inverse normal transformation to the original data
      fit_transformed = lm(formula = transformed ~ A) # perform the two sample t test again using the trasformed data
      p_transformed[i] = summary(fit_transformed)$coefficients[2,4] # collect the p values
      if (p_transformed[i] < siglevel) {
        error_transformed[i] = 1 # mark the errors in testing the transformed data
      }
    }
  }

  if (distribution == "exp") {
    rate = as.integer(readline(prompt = "rate : "))
    plot(density(rexp(n, rate)), main = "Density Curve")
    for(i in 1 : simulation){
      original = rexp(n, rate)
      fit_original = lm(formula = original ~ A)
      p_original[i] = summary(fit_original)$coefficients[2,4]
      if (p_original[i] < siglevel) {
        error_original[i] = 1
      }
      transformed = INT_function(original, n, k)
      fit_transformed = lm(formula = transformed ~ A)
      p_transformed[i] = summary(fit_transformed)$coefficients[2,4]
      if (p_transformed[i] < siglevel) {
        error_transformed[i] = 1
      }
    }
  }

  if (distribution == "chisq") {
    df = as.integer(readline(prompt = "df : "))
    plot(density(rchisq(n, df, ncp = 0)), main = "Density Curve")
    for(i in 1 : simulation){
      original = rchisq(n, df, ncp = 0)
      fit_original = lm(formula = original ~ A)
      p_original[i] = summary(fit_original)$coefficients[2,4]
      if (p_original[i] < siglevel) {
        error_original[i] = 1
      }
      transformed = INT_function(original, n, k)
      fit_transformed = lm(formula = transformed ~ A)
      p_transformed[i] = summary(fit_transformed)$coefficients[2,4]
      if (p_transformed[i] < siglevel) {
        error_transformed[i] = 1
      }
    }
  }

  if (distribution == "laplace") {
    location = as.integer(readline(prompt = "location : "))
    scale = as.integer(readline(prompt = "scale : "))
    plot(density(rlaplace(n, location, scale)), main = "Density Curve")
    for(i in 1 : simulation){
      original = rlaplace(n, location, scale)
      fit_original = lm(formula = original ~ A)
      p_original[i] = summary(fit_original)$coefficients[2,4]
      if (p_original[i] < siglevel) {
        error_original[i] = 1
      }
      transformed = INT_function(original, n, k)
      fit_transformed = lm(formula = transformed ~ A)
      p_transformed[i] = summary(fit_transformed)$coefficients[2,4]
      if (p_transformed[i] < siglevel) {
        error_transformed[i] = 1
      }
    }
  }

  if (distribution == "rayleigh") {
    scale = as.integer(readline(prompt = "scale : "))
    plot(density(rrayleigh(n, scale)), main = "Density Curve")
    for(i in 1 : simulation){
      original = rrayleigh(n, scale)
      fit_original = lm(formula = original ~ A)
      p_original[i] = summary(fit_original)$coefficients[2,4]
      if (p_original[i] < siglevel) {
        error_original[i] = 1
      }
      transformed = INT_function(original, n, k)
      fit_transformed = lm(formula = transformed ~ A)
      p_transformed[i] = summary(fit_transformed)$coefficients[2,4]
      if (p_transformed[i] < siglevel) {
        error_transformed[i] = 1
      }
    }
  }

  if (distribution == "weibull") {
    shape = as.integer(readline(prompt = "shape : "))
    scale = as.integer(readline(prompt = "scale : "))
    plot(density(rweibull(n, shape, scale)), main = "Density Curve")
    for(i in 1 : simulation){
      original = rweibull(n, shape, scale)
      fit_original = lm(formula = original ~ A)
      p_original[i] = summary(fit_original)$coefficients[2,4]
      if (p_original[i] < siglevel) {
        error_original[i] = 1
      }
      transformed = INT_function(original, n, k)
      fit_transformed = lm(formula = transformed ~ A)
      p_transformed[i] = summary(fit_transformed)$coefficients[2,4]
      if (p_transformed[i] < siglevel) {
        error_transformed[i] = 1
      }
    }
  }

  if (distribution == "cauchy") {
    location = as.integer(readline(prompt = "location : "))
    scale = as.integer(readline(prompt = "scale : "))
    plot(density(rcauchy(n, location, scale)), main = "Density Curve")
    for(i in 1 : simulation){
      original = rcauchy(n, location, scale)
      fit_original = lm(formula = original ~ A)
      p_original[i] = summary(fit_original)$coefficients[2,4]
      if (p_original[i] < siglevel) {
        error_original[i] = 1
      }
      transformed = INT_function(original, n, k)
      fit_transformed = lm(formula = transformed ~ A)
      p_transformed[i] = summary(fit_transformed)$coefficients[2,4]
      if (p_transformed[i] < siglevel) {
        error_transformed[i] = 1
      }
    }
  }

  if (distribution == "lognorm") {
    meanlog = as.integer(readline(prompt = "meanlog : "))
    sdlog = as.integer(readline(prompt = "sdlog : "))
    plot(density(rlnorm(n, meanlog, sdlog)), main = "Density Curve")
    for(i in 1 : simulation){
      original = rlnorm(n, meanlog, sdlog)
      fit_original = lm(formula = original ~ A)
      p_original[i] = summary(fit_original)$coefficients[2,4]
      if (p_original[i] < siglevel) {
        error_original[i] = 1
      }
      transformed = INT_function(original, n, k)
      fit_transformed = lm(formula = transformed ~ A)
      p_transformed[i] = summary(fit_transformed)$coefficients[2,4]
      if (p_transformed[i] < siglevel) {
        error_transformed[i] = 1
      }
    }
  }

  print(data.frame("TypeI Error Rate Original" = (length(which(error_original == 1)))/simulation, "TypeI Error Rate Transformed" = (length(which(error_transformed == 1)))/simulation)) # form the table with the type I error rates of the two tests using original and transformed data
  hist(p_original, main = "Distribution of p-values for Original Data")
  hist(p_transformed, main = "Distribution of p-values for Transformed Data")
}

