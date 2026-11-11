#' @title Modified Mann-Kendall Test For Serially Correlated Data Using the Yue and Wang (2004) Variance Correction Approach
#'
#' @description Time series data is often influenced by serial correlation. When data are not random and influenced by autocorrelation, modified Mann-Kendall tests may be used for trend detction. Yue and Wang (2004) have proposed variance correction approach to address the issue of serial correlation in trend analysis. Data are initially detrended and the effective sample size is calculated from the sample serial correlation coefficients of the detrended series at lags 1 to \code{lag.max}. No significance screening is applied to these coefficients.
#'
#' @importFrom stats acf median pnorm qnorm
#'
#' @param  x  - Time series data vector
#'
#' @param  lag.max  - Maximum lag of the serial correlation coefficients used in calculating the effective sample size. The default, \code{NULL}, uses all lags 1 to n-1 (the behaviour of earlier versions). See Details.
#'
#' @return  Corrected Zc  - Z statistic after variance Correction
#'
#' new P.value  - P-value after variance correction
#'
#' N/N*  - Effective sample size
#'
#' Original Z  - Original Mann-Kendall Z statistic
#'
#' Old P-value  - Original Mann-Kendall p-value
#'
#' Tau  - Mann-Kendall's Tau
#'
#' Sen's Slope  - Sen's slope
#'
#' old.variance - Old variance before variance Correction
#'
#' new.variance - Variance after correction
#'
#' @references Kendall, M. (1975). Rank Correlation Methods. Griffin, London, 202 pp.
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3): 245-259.
#'
#' @references Sen, P. K. (1968). Estimates of the Regression Coefficient Based on Kendall’s Tau. Journal of the American statistical Association, 63(324): 1379. <doi:10.2307/2285891>
#'
#' @references Yue, S. and Wang, C. Y. (2004). The Mann-Kendall test modified by effective sample size to detect trend in serially correlated hydrological series. Water Resources Management, 18(3): 201–218. <doi:10.1023/B:WARM.0000043140.61082.60>
#'
#' @references O'Brien, N. L., Burn, D. H., Annable, W. K. and Thompson, J. (2021). Water Resources Research, 57, e2020WR028886. <doi:10.1029/2020WR028886>
#'
#' @references Militino, A. F., Moradi, M. and Ugarte, M. D. (2020). Remote Sensing, 12(6), 1008. <doi:10.3390/rs12061008>
#'
#' @references Morrissey, M. L. (2026). Preprint, Zenodo. <doi:10.5281/zenodo.24015743>
#'
#' @details The variance correction approach suggested by Yue and Wang (2004) is implemeted in this function. By default (\code{lag.max = NULL}), sample serial correlation coefficients for all lags 1 to n-1 are used in calculating the effective sample size.
#'
#' Caution: with all lags and no significance screening, the correction factor n/n* is strongly biased low (for a series with no trend and no autocorrelation its expected value is about 1/3 before detrending, and smaller after Sen detrending), so the corrected test rejects the null hypothesis far too often. In simulations of trend-free white noise, about 45-50 percent of series are declared significant at the 5 percent level (O'Brien et al., 2021; Militino et al., 2020; Morrissey, 2026). With \code{lag.max = 1}, or with \code{mmky1lag}, the white-noise rejection rate is close to nominal (about 6-9 percent for n = 30-100); each additional lag inflates it further. If the estimated n/n* is not positive, the corrected variance is undefined and NaN is returned.
#'
#' @examples x<-c(Nile)
#' mmky(x)
#' mmky(x, lag.max = 1)
#'
#' @export
#'
mmky <-function(x, lag.max = NULL) {
  # Initialize the test parameters
  options(scipen = 999)
  # Time series Vector
  x = x
  # Modified Z statistic after variance correction as per Yue and Wang (2004) method
  z = NULL
  # Original Z statistic for Mann-Kendall test before variance correction
  z0 = NULL
  # Modified Z statistic after variance correction as per Yue and Wang (2004) method
  pval = NULL
  # Original p-value for Mann-Kendall test before variance correction
  pval0 = NULL
  # Initialize Mann-Kendall S statistic
  S = 0
  # Initialize Mann-Kendall Tau
  Tau = NULL
  # Correction factor n/n* value as per Yue and Wang (2004) method
  essf = NULL

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

  # Maximum lag used in the effective sample size (default: all lags 1 to n-1)

  if (is.null(lag.max)) {
    lag.max <- n - 1
  }
  if (length(lag.max) != 1 || !is.numeric(lag.max) || !is.finite(lag.max) ||
      lag.max < 1 || lag.max != round(lag.max)) {
    stop("lag.max must be NULL or a single positive integer")
  }
  if (lag.max > n - 1) {
    lag.max <- n - 1
    warning("lag.max is larger than n-1; n-1 lags are used")
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

  # Calculating trend-free series

  t=1:length(x)
  xn=(x[1:n])-((slp)*(t))

  # Calculating Mann-Kendall S statistic

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      S = S + sign(x[j]-x[i])
    }
  }

  # Calculating autocorrelation function of the detrended series (ro) for lags 1 to lag.max
  # Note: no significance screening is applied

  acf(xn, lag.max=lag.max, plot=FALSE)$acf[-1] -> ro

  # Calculating sum((1-(k/n))*ro[k]) for k=1,2,...,lag.max
  ess=0
  for(k in 1:lag.max)
      {
    ess=ess+(1-(k/n))*ro[k]
    }

  # Calculating variance correction factor (n/n*) as per Yue and Wang (2004)

  essf = 1 + 2*(ess)

  # Calculating Mann-Kendall variance before correction (Var(s))

  var.S = n*(n-1)*(2*n+5)*(1/18);
  if(length(unique(x)) < n) {
    unique(x) -> aux
    for (i in 1:length(aux)) {
      length(which(x == aux[i])) -> tie
      if (tie > 1) {
        var.S = var.S - tie*(tie-1)*(2*tie+5)*(1/18)
      }
    }
  }

  # Calculating new variance  Var(s)*=(Var(s))*(n/n*) as per Yue and Wang (2004)

  VS = var.S * essf

  # Calculating Z statistic values before and after variance correction

  if (S == 0) {
    z = 0
    z0 = 0
  }else
  if (S > 0) {
    z = (S-1)/sqrt(VS)
    z0 = (S-1)/sqrt(var.S)
  } else {
    z = (S+1)/sqrt(VS)
    z0 = (S+1)/sqrt(var.S)
  }

  # Calculating p-value before and after variance corection

  pval = 2*pnorm(-abs(z))
  pval0 = 2*pnorm(-abs(z0))

  # Calculating Kendall's Tau

  Tau = S/(.5*n*(n-1))

  # Listing all outputs

  return(c("Corrected Zc" = z,
           "new P-value" = pval,
           "N/N*" = essf,
           "Original Z" = z0,
           "old P.value" = pval0,
           "Tau" = Tau,
           "Sen's slope" = slp,
           "old.variance"=var.S,
           "new.variance"= VS))
}
