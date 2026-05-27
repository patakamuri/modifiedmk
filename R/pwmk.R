#' @title Mann-Kendall Test of Prewhitened Time Series Data in Presence of Serial Correlation Using the von Storch (1995) Approach
#'
#' @description When time series data are not random and influenced by autocorrelation, prewhitening the time series prior to application of trend test is suggested.
#'
#' @importFrom stats acf median pnorm qnorm
#'
#' @usage pwmk(x)
#'
#' @param  x  - Time series data vector
#'
#' @return  Z-Value  - Z statistic after prewhitening
#'
#' Sen's Slope  - Sen's slope for prewhitened series
#'
#' old. Sen's Slope  - Sen's slope for original data series (x)
#'
#' P-value  - P-value after prewhitening
#'
#' S  - Mann-Kendall S statistic
#'
#' Var(s) - Variance of S
#'
#' Tau  - Mann-Kendall's Tau
#'
#' @references Kendall, M. (1975). Rank Correlation Methods. Griffin, London, 202 pp.
#'
#' @references Kulkarni, A. and H. von Storch. 1995. Monte carlo experiments on the effects of serial correlation on the MannKendall test of trends. Meteorologische Zeitschrift N.F, 4(2): 82-85.
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3): 245-259.
#'
#' @references Salas, J.D. (1980). Applied modeling of hydrologic times series. Water Resources Publication, 484 pp.
#'
#' @references von Storch, V. H. (1995). Misuses of statistical analysis in climate research, In: Analysis of Climate Variability: Applications of Statistical Techniques, ed. von H. V. Storch and A. Navarra A. Springer-Verlag, Berlin: 11-26.
#'
#' @references Yue, S. and Wang, C. Y. (2002). Applicability of prewhitening to eliminate the influence of serial correlation on the Mann-Kendall test. Water Resources Research, 38(6), <doi:10.1029/2001WR000861>
#'
#' @details The lag-1 serial correlation coefficient is used for prewhitening.
#'
#' @examples x<-c(Nile)
#' pwmk(x)
#'
#' @export
#'
pwmk <-function(x) {
  # Initialize the test parameters
  options(scipen = 999)
  # Time series vector
  x = x
  # Modified Z statistic after prewhitening
  z = NULL
  # Modified p-value after prewhitening
  pval = NULL
  # Initialize Mann-Kendall S statistic
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

  n <- length(x)
  #Specify minimum input vector length
  if (n < 3) {
    stop("Input vector must contain at least three values")
  }

  # Calculating lag-1 autocorrelation coefficient (ro)

  acf(x, lag.max=1, plot=FALSE)$acf[-1] -> ro


  # Calculating prewhitened series

  a=1:(length(x)-1)
  b=2:(length(x))
  xn<-(x[b]-(x[a]*ro))

  n<-length(xn)
  n1<-length(x)


  # Calculating Mann-Kendall S statistic

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      S = S + sign(xn[j]-xn[i])
    }
  }

  # Calculating Mann-Kendall variance (Var(s))

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

  # Calculating Z statistic

  if (S == 0) {
    z = 0
  }else
  if (S > 0) {
    z = (S-1)/sqrt(var.S)
    } else {
    z = (S+1)/sqrt(var.S)
   }

  # Calculating p-value

  pval = 2*pnorm(-abs(z))


  # Calculating Kendall's Tau

  Tau = S/(.5*n*(n-1))


  # Calculating Sen's slope for original series

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


  return(c("Z-Value" = z,
           "Sen's Slope"= slp1,
           "old. Sen's Slope"= slp,
           "P-value" = pval,
           "S" = S,
           "Var(S)" = var.S,
           "Tau"=Tau))
}

