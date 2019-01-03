#' @title Modified Mann-Kendall Test For Serially Correlated Data Using the Hamed and Rao (1998) Variance Correction Approach Considering Only the First Three Lags
#'
#' @description Time series data is often influenced by serial correlation. When data are not random and influenced by autocorrelation, modified Mann-Kendall tests may be used for trend detction. Hamed and Rao (1998) have proposed variance correction approach to address the issue of serial correlation in Trend analysis. Data are initially detrended and the effective sample size is calulated using the ranks of significant serial correlation coefficients which are then used to correct the inflated (or deflated) variance of the test statistic.
#'
#' @importFrom stats acf median pnorm qnorm
#'
#' @usage mmkh3lag(x,ci=0.95)
#'
#' @param  x  - Time series data vector
#'
#' @param  ci - Confidence interval
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
#' @references Hamed, K. H. and Rao, A. R. (1998). A modified Mann-Kendall trend test for autocorrelated data. Journal of Hydrology, 204(1–4): 182–196. <doi:10.1016/S0022-1694(97)00125-X>.
#'
#' @references Kendall, M. (1975). Rank Correlation Methods. Griffin, London, 202 pp.
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3): 245-259.
#'
#' @references Rao, A. R., Hamed, K. H., & Chen, H.-L. (2003). Nonstationarities in hydrologic and environmental time series. Ringgold Inc., Portland, Oregon, 362 pp. <doi:10.1007/978-94-010-0117-5>
#'
#' @references Salas, J.D. (1980). Applied modeling of hydrologic times series. Water Resources Publication, 484 pp.
#'
#' @references Sen, P. K. (1968). Estimates of the Regression Coefficient Based on Kendall’s Tau. Journal of the American statistical Association, 63(324): 1379. <doi:10.2307/2285891>
#'
#' @details A detrended time series is constructed using Sen's slope and the lag-1 autocorreltation coefficient of the ranks of the data. The variance correction approach proposed by Hamed and Rao (1998) uses only significant lags of autocorrelation coefficients. As suggested by Rao et al. (2003), only the first three autocorrelation coefficients are used in this function.
#'
#' @examples x<-c(Nile)
#' mmkh3lag(x)
#'
#' @export

mmkh3lag <-function(x, ci=0.95) {
  # Initialize the test parameters

  # Time series vector
  x = x
  # Modified Z statistic after variance correction by Hamed Rao (1998) method
  z = NULL
  # Original Z statistic for Mann-Kendall test before variance correction
  z0 = NULL
  # Modified Z statistic after variance correction by Hamed and Rao (1998) method
  pval = NULL
  # Original p-value for Mann-Kendall test before variance correction
  pval0 = NULL
  # Initialize Mann-Kendall S statistic
  S = 0
  # Initialize Mann-Kendall Tau
  Tau = NULL
  # Correction factor n/n* value
  essf = NULL
  # Confidance interval
  ci = ci

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
  xn<-(x[1:n])-((slp)*(t))

  # Calculating Mann-Kendall S statistic

  n <- length(x)
    for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      S = S + sign(x[j]-x[i])
    }
  }

  # Calculating autocorrelation function of the ranks of observations (ro)

  acf(rank(xn), lag.max=3, plot=FALSE)$acf[-1] -> ro

  # Calculating significant autocorrelation at given confidance interval (rof)

  qnorm((1+ci)/2)/sqrt(n) -> sig
  rep(NA,length(ro)) -> rof
  for (i in 1:(length(ro))) {
    if(ro[i] > sig || ro[i] < -sig) {
      rof[i] <- ro[i]
    } else {
      rof[i] = 0
    }
  }

  # Calculating 2/(n*(n-1)*(n-2))

  2 / (n*(n-1)*(n-2)) -> cte

 # Calculating sum(((n-i)*(n-i-1)*(n-i-2)*rof[i]

  ess=0
      for (i in 1:3) {
            ess = ess + (n-i)*(n-i-1)*(n-i-2)*rof[i]
      }

  # Calculating variance correction factor (n/n*) as per Hamed and Rao (1998)

  essf = 1 + ess*cte

  # Calculating Mann-Kendall variance before correction (Var(s))

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

  # Calculating new variance  Var(s)*=(Var(s))*(n/n*) as per Hamed and Rao (1998)

  VS = var.S * essf

  # Calculating Z statistic values before and after variance correction

  if (S == 0) {
     z = 0
     z0 = 0
   }
   if (S > 0) {
     z = (S-1)/sqrt(VS)
     z0 = (S-1)/sqrt(var.S)
   } else {
     z = (S+1)/sqrt(VS)
     z0 = (S+1)/sqrt(var.S)
   }

   # Calculating p-value before and after variance correction

   pval = 2*pnorm(-abs(z))
   pval0 = 2*pnorm(-abs(z0))

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
