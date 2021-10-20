#' Inverse Normal Transformation
#'
#' @param x is the sample set.
#' @param n is the number of the samples.
#' @param k is an offset to ensure that all fractional ranks are strictly between zero and one.
#' @details Perform a rank-based inverse normal transformation.
#' @return The vector of the transformed data (Z-scores).

INT_function = function(x, n, k) {
  output = 0
  ranks = rank(x)
  for (i in 1:n) {
    output[i] = probitlink((ranks[i] - k)/(n + 1 - 2*k))
  }
  return(output)
}
