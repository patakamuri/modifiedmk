#' @title Mann-Kendall Trend Test of Time series Data with out any modifications
#'
#' @description Mann-Kendall trend test is a non-parametric trend test used to identify the monotonic trends present in time series data.
#'
#' @importFrom stats acf median pnorm qnorm
#'
#' @usage mkttest(x)
#'
#' @param  x  - Time series data vector
#'
#' @return  Z  - Mann- Kendall Z-statistic after variance Correction
#'
#' sen's slope  - sen's slope
#'
#' s  - Mann-Kendall 's'- statistic
#'
#' Var(s) - Variance of 's'
#'
#' P-value  - Mann-Kendall P-Value
#'
#' Tau  - Mann-Kendall's Tau
#'
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3), 245–259. <doi:10.1017/CBO9781107415324.004>
#'
#' @references Kendall, M. (1975). Multivariate analysis. Charles Griffin. Londres. 0-85264-234-2.
#'
#' @references sen, P. K. (1968). Estimates of the Regression Coefficient Based on Kendall’s Tau. Journal of the American statistical Association, 63(324), 1379. <doi:10.2307/2285891>
#'
#' @details Mann-Kendall trend test is a non-parametric trend tests which assuems no distribution of the data. Null hypothesis of the test is that there is no trend in the data and the alternate hypothesis is that the data represents a monotonic trend.
#'
#' @examples x<-c(Nile)
#' mkttest(x)
#'
#' @export
#'
mkttest <-function(x) {
  # Initialize the test Parameters

  # Time-Series Vector
  x = x
  # Mann-Kendall Z-Statistic
  z = NULL
  # Mann-Kendall P-value
  pval = NULL
  # Mann-Kendall 'S'- Statistic
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

  # Calculating Mann-Kendall 'S'- Statistic

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      S = S + sign(x[j]-x[i])
    }
  }

  # Calculating Mann-Kendall Variance (Var(s))

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

  # Calculating Z-Statistic values before and after Variance coorection

  if (S == 0) {
    z = 0
  }
  if (S > 0) {
    z = (S-1)/sqrt(var.S)
  } else {
    z = (S+1)/sqrt(var.S)
  }

  # Calculating P-Value before and after Variance coorection

  pval = 2*pnorm(-abs(z))


  # Calculating kendall's Tau

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


  return(c("Z-Value" = z,"Sen's slope"= slp, "P-value" = pval,"S" = S, "Var(S)" = var.S, "Tau"=Tau))

}

