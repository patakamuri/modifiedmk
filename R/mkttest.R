#' @title Mann-Kendall Trend Test of Time Series Data Without Modifications
#'
#' @description The Mann-Kendall trend test is a nonparametric trend test used to identify monotonic trends present in time series data.
#'
#' @importFrom stats acf median pnorm qnorm
#'
#' @usage mkttest(x)
#'
#' @param  x  - Time series data vector
#'
#' @return  Z  - Mann-Kendall Z statistic
#'
#' Sen's slope  - Sen's slope
#'
#' S  - Mann-Kendall S statistic
#'
#' Var(s) - Variance of S
#'
#' P-value  - Mann-Kendall p-value
#'
#' Tau  - Mann-Kendall's Tau
#'
#' @references Kendall, M. (1975). Rank Correlation Methods. Griffin, London, 202 pp.
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3): 245-259.
#'
#' @references Sen, P. K. (1968). Estimates of the Regression Coefficient Based on Kendall’s Tau. Journal of the American statistical Association, 63(324): 1379. <doi:10.2307/2285891>
#'
#' @details The Mann-Kendall trend test is a nonparametric trend tests which assumes no distribution of the data. The null hypothesis of the test is that there is no trend in the data and the alternative hypothesis is that the data represents a monotonic trend.
#'
#' @examples x<-c(Nile)
#' mkttest(x)
#'
#' @export
#'
mkttest <-function(x) {
  # Initialize the test parameters

  options(scipen = 999)

  # Time series vector
  x = x
  # Mann-Kendall Z statistic
  z = NULL
  # Mann-Kendall p-value
  pval = NULL
  # Mann-Kendall S statistic
  S = 0
  # Mann-Kendall var.S
  var.S = NULL
  # Mann-Kendall Tau
  Tau = NULL

  # To test whether the data is in vector format

  if (is.vector(x) == FALSE) {
    stop("Input data must be a vector")
  }

  # To test whether the data values are finite numbers and attempting to eliminate non-finite numbers

  if (any(is.finite(x) == FALSE)) {
    x[-c(which(is.finite(x) == FALSE))] -> x
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }

  n<-length(x)

  #Specify minimum input vector length
  if (n < 3) {
    stop("Input vector must contain at least three values")
  }

  # Calculating Mann-Kendall S statistic

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      S = S + sign(x[j]-x[i])
    }
  }

  # Calculating Mann-Kendall variance (Var(s))

  var.S = n*(n-1)*(2*n+5)*(1/18)
  if(length(unique(x)) < n) {
    unique(x) -> aux
    for (i in 1:length(aux)) {
      length(which(x == aux[i])) -> tie
      if (tie > 1) {
        var.S = var.S - tie*(tie-1)*(2*tie+5)*(1/18)
      }
    }
  }

  # Calculating Z statistic values

  if (S == 0) {
    z = 0
  } else
  if (S > 0) {
    z = (S-1)/sqrt(var.S)
  } else {
    z = (S+1)/sqrt(var.S)
  }

  # Calculating p-value

  pval = 2*pnorm(-abs(z))


  # Calculating Kendall's Tau

  Tau = S/(.5*n*(n-1))

  # Calculating Sen's slope

  rep(NA, n * (n - 1)/2) -> V
  k = 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      k = k+1
      V[k] = (x[j]-x[i])/(j-i)
    }
  }
  median(V,na.rm=TRUE)->slp


  return(c("Z-Value" = z,
           "Sen's slope"= slp,
           "S" = S,
           "Var(S)" = var.S,
           "P-value" = pval,
           "Tau"=Tau))

}

