#' @title Modified Mann-Kendall Test For Serially Correlated Data Using the Yue and Wang (2004) Variance Correction Approach
#'
#' @description Time series data is often influenced by serial correlation. When data are not random and influenced by autocorrelation, modified Mann-Kendall tests may be used for trend detction. Yue and Wang (2004) have proposed variance correction approach to address the issue of serial correlation in trend analysis. Data are initially detrended and the effective sample size is calculated using significant serial correlation coefficients.
#'
#' @importFrom stats acf median pnorm qnorm
#'
#' @param  x  - Time series data vector
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
#' @details The variance correction approach suggested by Yue and Wang (2004) is implemeted in this function. Serial correlation coefficients for all lags are used in calculating the effective sample size.
#'
#' @examples x<-c(Nile)
#' mmky(x)
#'
#' @export
#'
mmky <-function(x) {
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

  # Calculating autocorrelation function of the observations (ro)

  #lag.max can be edited to include greater numbers of lags

  acf(xn, lag.max=(n-1), plot=FALSE)$acf[-1] -> ro

  # Calculating significant autocorrelation at given confidance interval (rof)


  rep(NA,length(ro)) -> rof
  for (i in 1:(length(ro))) {
    rof[i] <- ro[i]
  }


  # Calculating sum(1-(k/n))*rof[i]) for k=1,2...,(n-1)
  ess=0
  for(k in 1:(n-1))
      {
    ess=ess+(1-(k/n))*rof[k]
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

