#' @title Nonparametric Block Bootstrapped Spearman's Rank Correlation Trend Test
#'
#' @description Significant serial correlation present in time series data can be accounted for using the nonparametric block bootstrap technique, which incorporates Spearman’s Rank Correlation trend test (Lehmann, 1975; Sneyers, 1990;Kundzewicz and Robson, 2000).  Predetermined block lengths are used in resampling the original time series, thus retaining the memory structure of the data.  If the value of the test statistic falls in the tails of the empirical bootstrapped distribution, there is likely a trend in the data.  The block bootstrap technique is powerful in the presence of autocorrelation (Khaliq et al. 2009; Önöz and Bayazit, 2012).
#'
#' @importFrom stats acf qnorm
#'
#' @importFrom boot tsboot boot.ci
#'
#' @usage bbssr(x, ci=0.95, nsim=2000, eta=1, bl.len=NULL)
#'
#' @param  x  - Time series data vector
#'
#' @param  ci - Confidence interval
#'
#' @param  nsim - Number of bootstrapped simulations
#'
#' @param  eta - Added to the block length
#'
#' @param bl.len - Block length
#'
#' @return Spearman's Correlation Coefficient - Spearman's correlation coefficient value
#'
#' Test Statistic - Z-transformed value to test significance \eqn{\rho(\sqrt{n-1})}
#'
#' Test Statistic Empirical Bootstrapped CI - Test statistic empirical bootstrapped confidence interval
#'
#' @references Box, G. E. P. and Jenkins, G. M. (1970). Time Series Analysis Forecasting and Control. Holden-Day, San Fransisco, California, 712 pp.
#'
#' @references Khaliq, M. N., Ouarda, T. B. M. J., Gachon, P., Sushama, L., and St-Hilaire, A. (2009). Identification of hydrological trends in the presence of serial and cross correlations: A review of selected methods and their application to annual flow regimes of Canadian rivers. Journal of Hydrology, 368: 117-130.
#'
#' @references Kundzewicz, Z. W. and Robson, A. J. (2000). Detecting Trend and Other Changes in Hydrological Data. World Climate Program-Water, Data and Monitoring. World Meteorological Organization, Geneva  (WMO/TD-No. 1013).
#'
#' @references Kundzewicz, Z. W. and Robson, A. J. (2004). Change detection in hydrological records-A review of the methodology. Hydrological Sciences Journal, 49(1): 7-19.
#'
#' @references Lehmann, E. L. (1975). Nonparametrics: statistical methods based on ranks. Holden-Day, Inc., California, 457 pp.
#'
#' @references Önöz, B. and Bayazit M. (2012). Block bootstrap for Mann-Kendall trend test of serially dependent data. Hydrological Processes, 26: 3552-3560.
#'
#' @references Sneyers, R. (1990). On the statistical analysis of series of observations. World Meteorological Organization, Technical Note no. 143, WMO no. 415, 192 pp.
#'
#' @references Svensson, C., Kundzewicz, Z. W., and Maurer, T. (2005). Trend detection in river flow series: 2. Floods and low-flow index series. Hydrological Sciences Journal, 50(5): 811-823.
#'
#' @details Block lengths are the number of contiguous significant serial correlations, to which the (\eqn{\eta}) term is added. A value of \eqn{\eta = 1} is used as the default as per Khaliq et al. (2009).  Alternatively, the user may define the block length.  2000 bootstrap replicates are recommended as per Svensson et al. (2005) and Önöz, B. and Bayazit (2012).
#'
#' @examples x<-c(Nile[1:10])
#' bbssr(x)
#'
#' @export
#'
bbssr<-function (x, ci = 0.95, nsim = 2000, eta = 1, bl.len = NULL)
{
  options(scipen = 999)
  x = x
  ci = ci
  nsim = nsim
  eta = eta
  bl.len = bl.len
  if (is.vector(x) == FALSE) {
    stop("Input data must be a vector")
  }
  n <- length(x)
  if (n < 4) {
    stop("Input vector must contain at least four values")
  }
  if (is.null(bl.len) == FALSE)
    if (bl.len > n) {
      stop("Block length must be less than the time series length")
    }
  if (any(is.finite(x) == FALSE)) {
    x <- x[-c(which(is.finite(x) == FALSE))]
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }
  n <- length(x)
  if (is.null(bl.len) == TRUE) {
    bd <- qnorm((1 + ci)/2)/sqrt(n)
    ro <- acf(x, lag.max = round(n/4), plot = FALSE)$acf[-1]
    sig.v <- rep(0, round(n/4))
    sig.vv <- rep(0, round(n/4))
    for (i in 1:round(n/4)) {
      if (-bd > ro[i] | bd < ro[i]) {
        sig.v[i] <- ro[i]
      }
    }
    if (all(sig.v == 0)) {
      min.sig <- 0
    }
    else {
      for (j in 1:length(sig.v)) {
        if (-bd > sig.v[j] | bd < sig.v[j]) {
          sig.vv[j]<-1
        }
      }
      min.sig.init <- rle(sig.vv)
      if (all(sig.vv == 0)) {
        min.sig <- 0
      } else {
        min.sig <- max(min.sig.init$lengths[min.sig.init$values != 0])
      }
    }
    bl.len <- min.sig + eta
  }
  SR.orig <- round(spear(x)[[1]], digits = 7)
  Z_trans <- round(spear(x)[[2]], digits = 7)
  SRfunc <- function(x) spear(x)[[2]]
  boot.out <- tsboot(x, SRfunc, R = nsim, l = bl.len, sim = "fixed")
  lb <- round(sort(boot.out$t)[(1 - ci)/2 * nsim], digits = 7)
  ub <- round(sort(boot.out$t)[(1 + ci)/2 * nsim], digits = 7)
  return(c("Spearman's Correlation Coefficient"=SR.orig,
           "Z-Transformed Test Statistic"=Z_trans,
           "Z-Transformed Test Statistic Empirical Bootstrapped CI Lower Bound"=lb,
           "Z-Transformed Test Statistic Empirical Bootstrapped CI Upper Bound"=ub))
}

