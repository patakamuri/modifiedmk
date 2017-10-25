#' @title Modified Mann-Kendall Test For Serially Correlated Data Using Yue and Wang (2004) Variance Correction Approach Considering Lag-1 Correlation Coefficient.
#'
#' @description Time series data is often influenced by serial-correlation. When data is not random and influenced by auto-correlation, Modified Mann-Kendall tests are to be used in trend detction. Yue and Wang (2004) have proposed variance correction approach to address the issue of serial correlation in Trend analysis. Trend is removed from the series and effective sample size is calculated using significant serial correlation coefficients.
#'
#' @importFrom stats acf median pnorm qnorm
#'
#' @param  x  - Time series data vector
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
#' @export
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3), 245–259. http://doi.org/10.1017/CBO9781107415324.004
#'
#' @references Kendall, M. (1975). Multivariate analysis. Charles Griffin. Londres. 0-85264-234-2.
#'
#' @references Sen, P. K. (1968). Estimates of the Regression Coefficient Based on Kendall’s Tau. Journal of the American Statistical Association, 63(324), 1379. http://doi.org/10.2307/2285891
#'
#' @references Yue, S., & Wang, C. Y. (2004). The Mann-Kendall test modified by effective sample size to detect trend in serially correlated hydrological series. Water Resources Management, 18(3), 201–218. http://doi.org/10.1023/B:WARM.0000043140.61082.60
#'
#' @examples x<-c(Nile)
#' mmky(x)
#'
#' @export
#'
mmky1lag <-function(x) {

# Initialize the test Parameters

# Time-Series Vector
  x = x
# Modified Z-Statistic after Variance Correction as per Yue and Wang (2004)) method
  z = NULL
# Original Z-Statistic for Mann-Kendall test before variance correction
  z0 = NULL
# Modified Z-Statistic after Variance Correction as per Yue and Wang (2004) method
  pval = NULL
# Original P-Value for Mann-Kendall test before variance correction
  pval0 = NULL
# Initialize Mann-Kendall 'S'- Statistic
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

  # Calculating Sen's slope

  n <- length(x)
  rep(NA, n * (n - 1)/2) -> V
  k = 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      k = k+1
      V[k] = (x[j]-x[i])/(j-i)
    }
  }
  median(V,na.rm=TRUE)->slp

# Calculating Trend-Free Series

  t=1:length(x)
  xn=(x[1:n])-((slp)*(t))

# Calculating Mann-Kendall 'S'- Statistic

    for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      S = S + sign(x[j]-x[i])
    }
  }

# Calculating auto-correlation function of the ranks of observations (ro)
  #lag.max can be edited to include more number of lags

  acf(xn, lag.max=1, plot=FALSE)$acf[-1] -> ro

# Calculating significant auto-correlation at given confidance interval (rof)


  rep(NA,length(ro)) -> rof
  for (i in 1:(length(ro))) {
    rof[i] <- ro[i]
  }


# Calculating sum(1-(k/n))*rof[i]) for k=1,2...,(n-1)

  ess=(1-(1/n))*(rof)


# Calculating variance correction factor (n/n*) as per Yue and Wang (2004)

  essf = 1 + 2*(ess)

# Calculating Mann-Kendall Variance before correction (Var(s))

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

# Calculating new variance  Var(s)*=(Var(s))*(n/n*)  as per Yue and Wang (2004)

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

# Listing all outputs



  return(c("Corrected Zc" = z, "new P-value" = pval,"Original Z" = z0, "old P.value" = pval0,"N/N*" = essf,"old.variance"=var.S, "new.variance"= VS))
}

