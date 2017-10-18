#' @title Modified Mann-Kendall Test For Serially Correlated Data Using Hamed and Rao (1998) Variance Correction Approach.
#'
#' @description Time series data is often influenced by serial-correlation. When data is not random and influenced by auto-correlation, Modified Mann-Kendall tests are to be used in trend detction. Hamed and Rao (1998) have proposed variance correction approach to address the issue of serial correlation in Trend analysis. Trend is removed from the series and effective sample size is calulated using significant serial correlation coefficients.
#'
#' @importFrom stats acf median pnorm qnorm
#'
#' @usage mmkh3lag(x,ci=0.95)
#'
#' @param  x  - Time series data vector
#'
#' @param  ci - Confidence Interval
#'
#' @return  Zc  - Z-Statistic after variance Correction
#'
#' new P.value  - P-Value after variance correction
#'
#' N/N*  - Effective sample size
#'
#' Z  - Original Mann-Kendall Z-Statistic
#'
#' P-value  - Original Mann-Kendall P-Value
#'
#' Tau  - Mann-Kendall's Tau
#'
#' Sen's Slope  - Sen's slope
#'
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3), 245–259. http://doi.org/10.1017/CBO9781107415324.004
#'
#' @references Kendall, M. (1975). Multivariate analysis. Charles Griffin. Londres. 0-85264-234-2.
#'
#' @references Sen, P. K. (1968). Estimates of the Regression Coefficient Based on Kendall’s Tau. Journal of the American Statistical Association, 63(324), 1379. http://doi.org/10.2307/2285891
#'
#' @references Hamed, K. H., & Ramachandra Rao, A. (1998). A modified Mann-Kendall trend test for autocorrelated data. Journal of Hydrology, 204(1–4), 182–196. http://doi.org/10.1016/S0022-1694(97)00125-X
#'
#' @references Rao, A. R., Hamed, K. H., & Chen, H.-L. (2003). Nonstationarities in hydrologic and environmental time series. http://doi.org/10.1007/978-94-010-0117-5
#'
#' @references Salas, J.D., 1980. Applied modeling of hydrologic times series. Water Resources Publication.
#'
#' @details Trend free time series is constructed by calculating Sen's slope and Auto Correlation coefficient AR(1). Variance correction approach proposed by Hamed and Rao (1998) uses only significant values from all the available values of Auto-Correlation Coefficients. As suggested by Rao, A. R., Hamed, K. H., & Chen, H.-L. (2003), only first three Auto-Correlation coefficients are used.
#'
#' @examples x<-c(Nile)
#' mmkh3lag(x)
#'
#' @export

mmkh3lag <-function(x, ci=0.95) {

# Initialize the test Parameters

# Time-Series Vector
  x = x
# Modified Z-Statistic after Variance Correction by Hamed&Rao(1998) method
  z = NULL
# Original Z-Statistic for Mann-Kendall test before variance correction
  z0 = NULL
# Modified Z-Statistic after Variance Correction by Hamed&Rao(1998) method
  pval = NULL
# Original P-Value for Mann-Kendall test before variance correction
  pval0 = NULL
# Initialize Mann-Kendall 'S'- Statistic
  S = 0
# Initialize Mann-Kendall Tau
  Tau = NULL
# Correction factor n/n* value
  essf = NULL
# Confidance Interval
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

# Calculating Mann-Kendall 'S'- Statistic

  n <- length(x)
    for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      S = S + sign(x[j]-x[i])
    }
  }

# Calculating auto-correlation function of the ranks of observations (ro)

  acf(rank(x), lag.max=3, plot=FALSE)$acf[-1] -> ro

# Calculating significant auto-correlation at given confidance interval (rof)

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

# Calculating variance correction factor (n/n*) as per Hamed and Rao(1994)

  essf = 1 + ess*cte

# Calculating Mann-Kendall Variance before correction (Var(s))

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

# Calculating new variance  Var(s)*=(Var(s))*(n/n*)  as per Hamed and Rao(1994)

  VS = var.S * essf

# Calculating Z-Statistic values before and after Variance coorection

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

# Calculating P-Value before and after Variance coorection

  pval = 2*pnorm(-abs(z))
      pval0 = 2*pnorm(-abs(z0))

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

  return(list("Corrected Zc" = z, "Corrected p.value" = pval,"N/N*s" = essf,"Z" = z0, "p.value" = pval0,  "tau" = Tau,  "Sen's Slope" = slp))
}

