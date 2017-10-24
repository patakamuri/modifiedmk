#' @title Mann-Kendall Test of Pre-Whitened Time Series Data in Presence of Serial Correlation Using (Yue and Wang 2002) Approach.
#'
#' @description When the time series data is not random and influenced by auto-correlation, Pre-Whitening the time series prior to application of trend test is suggested.
#'
#' @importFrom stats acf median pnorm qnorm
#'
#' @usage pwmk(x)
#'
#' @param  x  - Time series data vector
#'
#' @return  Z-Value  - Z-Statistic after variance Correction
#'
#' Sen's Slope  - Sen's slope for Prewhitened series
#'
#' old. Sen's Slope  - Sen's slope for Original data series 'x'
#'
#' P-value  - P-Value after variance correction
#'
#' S  - Mann-Kendall 'S'- statistic
#'
#' Var(s) - Variance of 's'
#'
#' Tau  - Mann-Kendall's Tau
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3), 245–259. <doi:10.1017/CBO9781107415324.004>
#'
#' @references Kendall, M. (1975). Multivariate analysis. Charles Griffin. Londres. 0-85264-234-2.
#'
#' @references Sen, P. K. (1968). Estimates of the Regression Coefficient Based on Kendall’s Tau. Journal of the American Statistical Association, 63(324), 1379. <doi:10.2307/2285891>
#'
#' @references Yue, S., & Wang, C. Y. (2002). Applicability of prewhitening to eliminate the influence of serial correlation on the Mann-Kendall test. Water Resources Research, 38(6), 4-1-4–7. <doi:10.1029/2001WR000861>
#'
#' @references Salas, J.D., (1980). Applied modeling of hydrologic times series. Water Resources Publication.
#'
#' @details Pre-Whitening involves calculating lag-1 serial correlation coefficient and calculating new-series.
#'
#' @examples x<-c(Nile)
#' pwmk(x)
#'
#' @export
#'
pwmk <-function(x) {

  # Initialize the test Parameters

  # Time-Series Vector
  x = x
  # Modified Z-Statistic after Pre-Whitening
  z = NULL
  # Modified P-value after Pre-Whitening
  pval = NULL
  # Initialize Mann-Kendall 'S'- Statistic
  S = 0
  # Initialize Mann-Kendall var.S
  var.S = NULL
  # Initialize Mann-Kendall Tau
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

  # Calculating lag-1 auto-correlation coefficient (ro)

  acf(x, lag.max=1, plot=FALSE)$acf[-1] -> ro


  # Calculating pre-whitened Series

  a=1:(length(x)-1)
  b=2:(length(x))
  xn<-(x[b]-(x[a]*ro))

  n<-length(xn)
  n1<-length(x)


  # Calculating Mann-Kendall 'S'- Statistic

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      S = S + sign(xn[j]-xn[i])
    }
  }

  # Calculating Mann-Kendall Variance (Var(s))

  var.S = n*(n-1)*(2*n+5)*(1/18)
  if(length(unique(xn)) < n) {
    unique(xn) -> aux
    for (i in 1:length(aux)) {
      length(which(xn == aux[i])) -> tie
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


  # Calculating Sen's slope for original series 'x'

  rep(NA, n1 * (n1 - 1)/2) -> V
  k = 0
  for (i in 1:(n1-1)) {
    for (j in (i+1):n1) {
      k = k+1
      V[k] = (x[j]-x[i])/(j-i)
    }
  }
  median(V,na.rm=TRUE)->slp

  # Calculating Sen's slope for PW series

  rep(NA, n * (n - 1)/2) -> W
  m = 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      m = m+1
      W[m] = (xn[j]-xn[i])/(j-i)
    }
  }
  median(W,na.rm=TRUE)->slp1


  return(list("Z-Value" = z,"Sen's Slope"= slp1, "old. Sen's Slope"= slp,"P.value" = pval,"S" = S, "Var(S)" = var.S, "Tau"=Tau))
}

