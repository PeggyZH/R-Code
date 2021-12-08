#' Power Analysis
#' @param distribution contains 4 types of distributions which are Normal Distribution, Chi-square Distribution, LaPlace Distribution, and Weibull Distribution. The particular distribution for which the data is generated from can be specified using this parameter.
#' @param simulation is the number of simulations.
#' @param n is the number of the samples.
#' @param k is an offset to ensure that all fractional ranks are strictly between zero and one.
#' @details Conduct a power analysis on the two sample t tests and the permutation tests performed on the data generated from the desired distribution and the same set of data after inverse normal transformation. Observe whether the inverse normal transformation is significant and whether the choice of test is influencial. (Reproducing Table 4 in Article "Rank-Based Inverse Normal Transformations are Increasingly Used, But are They Merited?")
#' @details Need to install INT_Function package first.
#' @return a table that contains the powers of the two sample t tests and the permutation tests under significance levels of 0.05, 0.01, 0.001 using the original data and the transformed data (inverse normal transformation) respectively.


pwanalysis <- function(simulation, n, k, distribution) {
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

  for (i in 1:simulation) { # number of simulations
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
  }

  # the above process is repeated but perform the permutation test rather than the two sample t test
  for (i in 1:simulation) {
    if (distribution == "norm") {
      x = rnorm(n, mean = 0, sd = 1)
      y = rnorm(n, mean = 0.5, sd = 1)
    }
    if (distribution == "laplace") {
      x = rlaplace(n, location = 0, scale = sqrt(0.5))
      y = rlaplace(n, location = 0.5, scale = sqrt(0.5))
    }
    if (distribution == "chisq") {
      x = rchisq(n, df = 1)
      y = rchisq(n, df = 1) + 0.5
    }
    if (distribution == "weibull") {
      x = rweibull(n, shape = 1, scale = 0.5)
      y = rweibull(n, shape = 1, scale = 0.5) + 0.5
    }
    x = rnorm(n, mean = 0, sd = 1)
    y = rnorm(n, mean = 0.5, sd = 1)
    ptest_original[i] = permTS(x, y)$p.value
    if (ptest_original[i] < 0.05) {
      ptest_original_005[i] = 1
    }
    if (ptest_original[i] < 0.01) {
      ptest_original_001[i] = 1
    }
    if (ptest_original[i] < 0.001) {
      ptest_original_0001[i] = 1
    }
    total = c(x, y)
    transformed = INT_function(total, 2*n, k)
    splited = split(transformed, ceiling(seq_along(transformed) / n))
    newx = splited$`1`
    newy = splited$`2`
    ptest_transformed[i] = permTS(newx, newy)$p.value
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
