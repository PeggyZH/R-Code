#' Inverse Normal Transformation
#'
#' @param x is the sample set.
#' @param n is the number of the samples.
#' @param k is an offset to ensure that all fractional ranks are strictly between zero and one.
#'
#' @import VGAM
#'
#' @return The vector of the transformed data (Z-scores).
#'
#' @details Perform a rank-based inverse normal transformation.
#'
#' @export

INT_function = function(x, n, k) {
  output = 0
  ranks = rank(x) # return the sample ranks of the values
  for (i in 1 : n) {
    output[i] = probitlink((ranks[i] - k)/(n + 1 - 2*k))
  }
  return(output)
}

#' Type I Error Rate and Distribution of P Values
#'
#' @param distribution contains 8 types of distributions which are Normal Distribution, Chi-square Distribution, Exponential Distribution, LaPlace Distribution, Rayleigh Distribution, Weibull Distribution, Cauchy Distribution, and Log Normal Distribution. The particular distribution for which the data is generated from can be specified using this parameter.
#' @param simulation is the number of simulations.
#' @param n is the number of the samples.
#' @param k is an offset to ensure that all fractional ranks are strictly between zero and one.
#' @param siglevel is the significance level for the two sample test.
#' @param param is the parameter for the distributions. Normal Distribution: param = c(mean, sd). Chi-square Distribution: param = degree of freedom. Exponential Distribution: param = rate. LaPlace Distribution: param = c(the location parameter which is the mean, the scale parameter which must consist of positive values). Rayleigh Distribution: prarm = the scale parameter. Weibull Distribution: param = c(shape, scale). Cauchy Distribution: param = c(location, scale). Log Normal Distribution: param = c(mean on the log scale, sd on the log scale).
#'
#' @details Use a given matrix A = [0, 0, 0 ...... 1, 1, 1] to fit the data generated from a target distribution. Then transform the data using the inverse normal transformation and fit it with the matrix A again. Calculate the type I error rates for the two tests and show the distributions of the p values to evaluate the significance of the transformation.
#'
#' @return a table of type I error rates of the two two sample tests conducted using original and transformed data respectively and two histograms that display the distributions of p values in the two tests.
#'
#' @export

dpvalue <- function(simulation, n, k, siglevel, distribution, param) {
X1 = c(rep(1, n/2), rep(0, n/2))
X2 = c(rep(0, n/2), rep(1, n/2))
A = matrix(c(X1, X2), nrow = n, ncol=2) # construct the matrix A = [0, 0, 0 ...... 1, 1, 1]
p_original = vector(length = simulation) # store the p values for original data
error_original = NA # mark the position where type I error occurs in original data
p_transformed = vector(length = simulation) # store the p values for transformed data
error_transformed = NA # mark the position where type I error occurs in transformed data

for(i in 1 : simulation){
  if (distribution == "norm") {
    mean = param[1] # input the mean of the distribution
    sd = param[2] # input the standard deviation of the distribution
    original = rnorm(n, mean, sd) # generate data from the desired distribution
  }
  if (distribution == "exp") {
    rate = param # input rate
    original = rexp(n, rate)
  }
  if (distribution == "chisq") {
    df = param # input degree of freedom
    original = rchisq(n, df, ncp = 0)
  }
  if (distribution == "laplace") {
    location = param[1] # input location
    scale = param[2] # input scale
    original = rlaplace(n, location, scale)
  }
  if (distribution == "rayleigh") {
    scale = param # input sclae
    original = rrayleigh(n, scale)
  }
  if (distribution == "weibull") {
    shape = param[1] # input shape
    scale = param[2] # input scale
    original = rweibull(n, shape, scale)
  }
  if (distribution == "cauchy") {
    location = param[1] # input location
    scale = param[2] # input scale
    original = rcauchy(n, location, scale)
  }
  if (distribution == "lognorm") {
    meanlog = param[1] # input mean on the log scale
    sdlog = param[2] # input standard deviation on the log scale
    original = rlnorm(n, meanlog, sdlog)
  }

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

plot(density(original), main = "Density Curve") # plot the density curve for reference
print(data.frame("TypeI Error Rate Original" = (length(which(error_original == 1)))/simulation, "TypeI Error Rate Transformed" = (length(which(error_transformed == 1)))/simulation)) # form the table with the type I error rates of the two tests using original and transformed data
hist(p_original, main = "Distribution of p-values for Original Data") # draw a histogram to display the distribution of the p values calculated from original data
hist(p_transformed, main = "Distribution of p-values for Transformed Data") # draw a histogram to display the distribution of the p values calculated from transformed data
}

#' Power Analysis
#'
#' @param distribution contains 4 types of distributions which are Normal Distribution, Chi-square Distribution, LaPlace Distribution, and Weibull Distribution. The particular distribution for which the data is generated from can be specified using this parameter.
#' @param simulation is the number of simulations.
#' @param n is the number of the samples.
#' @param k is an offset to ensure that all fractional ranks are strictly between zero and one.
#'
#' @import perm
#'
#' @details Conduct a power analysis on the two sample t tests and the permutation tests performed on the data generated from the desired distribution and the same set of data after inverse normal transformation. Observe whether the inverse normal transformation is significant and whether the choice of test is influencial. (Reproducing Table 4 in Article "Rank-Based Inverse Normal Transformations are Increasingly Used, But are They Merited?")
#'
#' @return a table that contains the powers of the two sample t tests and the permutation tests under significance levels of 0.05, 0.01, 0.001 using the original data and the transformed data (inverse normal transformation) respectively.
#'
#' @export

pwanaly <- function(simulation, n, k, distribution) {
  ttest_original = vector(length = simulation) # store the p values for t tests using original data

  # store the position where the p value is less than the significance level for t tests using original data
  ttest_original_005 = NA # siglevel = 0.05
  ttest_original_001 = NA # siglevel = 0.01
  ttest_original_0001 = NA # siglevel = 0.001

  ttest_transformed = vector(length = simulation) # store the p values for t tests using transformed data

  # store the position where the p value is less than the significance level for t tests using transformed data
  ttest_transformed_005 = NA # siglevel = 0.05
  ttest_transformed_001 = NA # siglevel = 0.01
  ttest_transformed_0001 = NA # siglevel = 0.001

  ptest_original = vector(length = simulation) # store the p values for permutation tests using original data

  # store the position where the p value is less than the significance level for permutation tests using original data
  ptest_original_005 = NA # siglevel = 0.05
  ptest_original_001 = NA # siglevel = 0.01
  ptest_original_0001 = NA # siglevel = 0.001

  ptest_transformed = vector(length = simulation) # store the p values for permutation tests using transformed data

  # store the position where the p value is less than the significance level for permutation tests using transformed data
  ptest_transformed_005 = NA # siglevel = 0.05
  ptest_transformed_001 = NA # siglevel = 0.01
  ptest_transformed_0001 = NA # siglevel = 0.001

  for (i in 1:simulation) {
    # construct two groups samples from the desired distribution
    # there is a 0.5 difference in group means
    if (distribution == "norm") {
      a = rnorm(n, mean = 0, sd = 1)
      b = rnorm(n, mean = 0.5, sd = 1)
    }
    if (distribution == "laplace") {
      a = rlaplace(n, location = 0, scale = sqrt(0.5))
      b = rlaplace(n, location = 0.5, scale = sqrt(0.5))
    }
    if (distribution == "chisq") {
      a = rchisq(n, df = 1)
      b = rchisq(n, df = 1) + 0.5
    }
    if (distribution == "weibull") {
      a = rweibull(n, shape = 1, scale = 0.5)
      b = rweibull(n, shape = 1, scale = 0.5) + 0.5
    }

    ttest_original[i] = t.test(a, b)$p.value # perform the two sample t test and collect the p values

    # mark the position where the p value is less than the significance level
    if (ttest_original[i] < 0.05) {
      ttest_original_005[i] = 1
    }
    if (ttest_original[i] < 0.01) {
      ttest_original_001[i] = 1
    }
    if (ttest_original[i] < 0.001) {
      ttest_original_0001[i] = 1
    }

    total = c(a, b) # combine the two groups into one
    transformed = INT_function(total, 2*n, k) # perform the inverse normal transformation
    splited = split(transformed, ceiling(seq_along(transformed) / n))
    newa = splited$`1`
    newb = splited$`2` # split the transformed data back to two groups of samples along the sequence

    ttest_transformed[i] = t.test(newa, newb)$p.value # perform the two sample t test and collect the p values

    # mark the position where the p value is less than the significance level
    if (ttest_transformed[i] < 0.05) {
      ttest_transformed_005[i] = 1
    }
    if (ttest_transformed[i] < 0.01) {
      ttest_transformed_001[i] = 1
    }
    if (ttest_transformed[i] < 0.001) {
      ttest_transformed_0001[i] = 1
    }

    # the above process is repeated but perform the permutation test rather than the two sample t test
    ptest_original[i] = permTS(a, b)$p.value
    if (ptest_original[i] < 0.05) {
      ptest_original_005[i] = 1
    }
    if (ptest_original[i] < 0.01) {
      ptest_original_001[i] = 1
    }
    if (ptest_original[i] < 0.001) {
      ptest_original_0001[i] = 1
    }
    ptest_transformed[i] = permTS(newa, newb)$p.value
    if (ptest_transformed[i] < 0.05) {
      ptest_transformed_005[i] = 1
    }
    if (ptest_transformed[i] < 0.01) {
      ptest_transformed_001[i] = 1
    }
    if (ptest_transformed[i] < 0.001) {
      ptest_transformed_0001[i] = 1
    }
  }

  # calculate the powers and form them into one data frame
  data = (data.frame(c("Power Original 0.05" = (length(which(ttest_original_005 == 1)))/simulation, "Power Original 0.01" = (length(which(ttest_original_001 == 1)))/simulation, "Power Original 0.001" = (length(which(ttest_original_0001 == 1)))/simulation, "Power Transformed 0.05" = (length(which(ttest_transformed_005 == 1)))/simulation, "Power Transformed 0.01" = (length(which(ttest_transformed_001 == 1)))/simulation, "Power Transformed 0.001" = (length(which(ttest_transformed_0001 == 1)))/simulation), c("Power Original 0.05" = (length(which(ptest_original_005 == 1)))/simulation, "Power Original 0.01" = (length(which(ptest_original_001 == 1)))/simulation, "Power Original 0.001" = (length(which(ptest_original_0001 == 1)))/simulation, "Power Transformed 0.05" = (length(which(ptest_transformed_005 == 1)))/simulation, "Power Transformed 0.01" = (length(which(ptest_transformed_001 == 1)))/simulation, "Power Transformed 0.001" = (length(which(ptest_transformed_0001 == 1)))/simulation)))
  # rename the columns and categorize the results into the ones using t test and the ones using permutation test
  colnames(data) = c("T Test", "Permutation Test")
  # display the result
  print(data)
}
