#' @title Mann-Kendall Trend Test Applied to Trend-Free Prewhitened Time Series Data in Presence of Serial Correlation Using Yue et al. (2002) Approach
#'
#' @description When the time series data are not random and influenced by autocorrelation, the trend component is removed from the data and is prewhitened prior to the application of the trend test.
#'
#' @importFrom stats acf median pnorm qnorm
#'
#' @usage tfpwmk(x)
#'
#' @param  x  - Time series data vector
#'
#' @return  Z-Value  - Z statistic after trend-free prewhitening (TFPW)
#'
#' Sen's Slope  - Sen's slope for TFPW series
#'
#' Old Sen's Slope  - Sen's slope for original data series (x)
#'
#' P-value  - P-value after trend-free prewhitening
#'
#' S  - Mann-Kendall S statistic
#'
#' Var(s) - Variance of S
#'
#' Tau  - Mann-Kendall's Tau
#'
#'
#' @references Kendall, M. (1975). Rank Correlation Methods. Griffin, London, 202 pp.
#'
#' @references Kulkarni, A. and H. von Storch. 1995. Monte carlo experiments on the effects of serial correlation on the MannKendall test of trends. Meteorologische Zeitschrift N.F, 4(2): 82-85.
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3): 245-259.
#'
#' @references Salas, J.D. (1980). Applied modeling of hydrologic times series. Water Resources Publication, 484 pp.
#'
#' @references Sen, P. K. (1968). Estimates of the Regression Coefficient Based on Kendall’s Tau. Journal of the American Statistical Association, 63(324): 1379. <doi:10.2307/2285891>
#'
#' @references von Storch, V. H. (1995). Misuses of statistical analysis in climate research, In: Analysis of Climate Variability: Applications of Statistical Techniques, ed. von H. V. Storch and A. Navarra A. Springer-Verlag, Berlin: 11-26.
#'
#' @references Yue, S., Pilon, P., Phinney, B., and Cavadias, G. (2002). The influence of autocorrelation on the ability to detect trend in hydrological series. Hydrological Processes, 16(9): 1807–1829. <doi:10.1002/hyp.1095>
#'
#' @details The linear trend component is removed from the original data and then prewhitened using the lag-1 serial correlation coefficient. The prewhitening data are then tested with Mann-Kendall trend test.
#'
#' @examples x<-c(Nile)
#' tfpwmk(x)
#'
#' @export
#'
tfpwmk <-function(x) {
  # Initialize the test parameters

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

  n<-length(x)
  
  #Specify minimum input vector length
  if (n < 3) {
    stop("Input vector must contain at least three values")
  }
  
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

  # Calculating trend-free series (xt)

  t=1:length(x)
  xt<-(x[1:n])-((slp)*(t))

  # Calculating lag-1 autocorrelation coefficient of trend-free series (ro)

  acf(xt, lag.max=1, plot=FALSE)$acf[-1] -> ro


  # Calculating Trend-Free Pre-Whitened Series (xp)

  a=1:(length(xt)-1)
  b=2:(length(xt))
  xp<-(xt[b]-(xt[a]*ro))

  # Calculating blended series to which MK test is to be applied in presence of significant serial correlation

  l<-length(xp)
  q=1:l
  y<-(xp[1:l]+((slp)*(q)))

  n1<-length(y)

  # Calculating Mann-Kendall S statistic

  for (i in 1:(n1-1)) {
    for (j in (i+1):n1) {
      S = S + sign(y[j]-y[i])
    }
  }

  # Calculating Mann-Kendall variance (Var(s))

  var.S = n1*(n1-1)*(2*n1+5)*(1/18)
  if(length(unique(y)) < n1) {
    unique(y) -> aux
    for (i in 1:length(aux)) {
      length(which(y == aux[i])) -> tie
      if (tie > 1) {
        var.S = var.S - tie*(tie-1)*(2*tie+5)*(1/18)
      }
    }
  }

  # Calculating Z statistic

  if (S == 0) {
    z = 0
  }
  if (S > 0) {
    z = (S-1)/sqrt(var.S)
  } else {
    z = (S+1)/sqrt(var.S)
  }

  # Calculating p-value

  pval = 2*pnorm(-abs(z))


  # Calculating Kendall's Tau

  Tau = S/(.5*n1*(n1-1))

  # Calculating Sen's slope for TFPW series

  rep(NA, n1 * (n1 - 1)/2) -> W
  m = 0
  for (i in 1:(n1-1)) {
    for (j in (i+1):n1) {
      m = m+1
      W[m] = (y[j]-y[i])/(j-i)
    }
  }
  median(W,na.rm=TRUE)->slp1

  return(c("Z-Value" = z,
           "Sen's Slope"= slp1,
           "Old Sen's Slope"= slp,
           "P-value" = pval,
           "S" = S,
           "Var(S)" = var.S,
           "Tau"=Tau))
}
