#' @title Block bootstrapped Mann-Kendall Trend Test
#'
#' @description Block lengths are automatically calculated using the lag of the smallest signigicant auto-correlation coefficient.  And eta value of 1 is used as a default value following the results of Khaliq et al. (2009), 2000 bootstrapped replicates are recommended following the results of Svensson et al. (2005).
#'
#' @importFrom stats acf qnorm
#'
#' @importFrom boot tsboot boot.ci
#'
#' @usage bbsmk(x, ci=0.95, nsim=2000, eta=1, bl.len=NULL)
#'
#' @param  x  - Time series data vector
#'
#' @param  ci - Confidence Interval
#'
#' @param  nsim - Number of simulations
#'
#' @param  eta - Added to the block length
#'
#' @param bl.len - Block length
#'
#' @return  Z  - Mann- Kendall Z-statistic
#'
#' slp  - sen's slope
#'
#' S  - Mann-Kendall 's'- statistic
#'
#' Tau  - Mann-Kendall's Tau
#'
#' lb- Lower bound of confidence interval value
#'
#' ub- Upper bound of confidence interval value
#'
#' @references Box, G. E. P. and Jenkins, G. M. (1970). Time Series Analysis Forecasting and Control. Holden-Day, San Fransisco, California, 712 pp.
#'
#' @references Kendall, M. (1975). Multivariate analysis. Charles Griffin. Londres. 0-85264-234-2.
#'
#' @references Khaliq, M. N., Ouarda, T. B. M. J., Gachon, P., Sushama, L., and St-Hilaire, A. (2009). Identification of hydrological trends in the presence of serial and cross correlations: A review of selected methods and their application to annual flow regimes of Canadian rivers. Journal of Hydrology, 368: 117-130.
#'
#' @references Kundzewicz, Z. W. and Robson, A. J. (2000). Detecting Trend and Other Changes in Hydrological Data. World Climate Program-Data and Monitoring. World Meteorological Organization, Geneva  (WMO/TD-No. 1013).
#'
#' @references Kundzewicz, Z. W. and Robson, A. J. (2004). Change detection in hydrological records-A review of the methodology. Hydrological Sciences Journal, 49(1): 7-19.
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3), 245-259. <doi:10.1017/CBO9781107415324.004>
#'
#' @references Önöz , B. and Bayazit M. (2012). Block bootstrap for Mann-Kendall trend test of serially dependent data. Hydrological Processes, 26: 3552-3560.
#'
#' @references Svensson, C., Kundzewicz, Z. W., and Maurer, T. (2005). Trend detection in river flow series: 2. Floods and low-flow index series. Hydrological Sciences Journal, 50(5): 811-823.
#'
#' @details The block bootstrap is used along with the non-parametric Mann-Kendall trend test.  A test statistic falling in the tails of the simulated empirical distribution, the results is likely significant.
#'
#' @examples x<-c(Nile)
#' bbsmk(x)
#'
#' @export
#'
bbsmk <- function(x,ci=0.95,nsim=2000,eta=1, bl.len=NULL) {
  # Initialize the test Parameters

  # Time-Series Vector
  x = x
  # Confidance Interval
  ci = ci
  #Number of Simulations
  nsim=nsim
  #Value of eta
  eta=eta
  #Initialization of block length
  bl.len=bl.len
  # Mann-Kendall Tau
  Tau = NULL
  # To test whether the data is in vector format

  if (is.vector(x) == FALSE) {
    stop("Input data must be a vector")
  }

  n<-length(x)

  #Specify minimum input vector length
  if (n < 4) {
    stop("Input vector must contain at least four values")
  }

  #Specify minimum block length
  if (is.null(bl.len) == FALSE)
    if (bl.len > n) {
      stop("Block length must be less than the time series length")
    }

  # To test whether the data values are finite numbers and attempting to eliminate non-finite numbers
  if (any(is.finite(x) == FALSE)) {
    x[-c(which(is.finite(x) == FALSE))] -> x
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }

  if (is.null(bl.len) == TRUE) {
    #bounds of the confidence intervals of the acf function

    bd <- qnorm((1 + ci)/2)/sqrt(n)

    # Calculating auto-correlation function of the observations (ro)

    ro <- acf(x, lag.max=round(n/4), plot=FALSE)$acf[-1]

    #Initialize vector of significant lags of auto-correlation

    sig.v <- rep(0, round(n/4))

    #Identify the number of contiguous significant serial correlations

    for (i in 1:round(n/4)) {
      if (-bd > ro[i] | bd < ro[i]) {
        sig.v[i]<-ro[i]
      }
    }

    if (all(sig.v == 0)) {
      min.sig<-0
    } else {
      min.sig.init<-rle(sig.v)
      min.sig<-max(min.sig.init$lengths[min.sig.init$values != 0])
    }

    #Block length

    bl.len <- min.sig + eta

  }

    #Block bootstrap using Mann Kendall

  MK.orig <- mkttest(x)
  Z<-round(MK.orig["Z-Value"], digits = 7)
  slp<-round(MK.orig["Sen's slope"], digits = 7)
  S<-MK.orig["S"]
  Tau <- MK.orig["Tau"]
  MKtau <- function(x) mkttest(x)[["Tau"]]
  boot.out.MKtau <- tsboot(x, MKtau, R=nsim, l=bl.len, sim="fixed")
  MKZ <- function(x) mkttest(x)[[1]]
  boot.out.Zval <- tsboot(x, MKZ, R=nsim, l=bl.len, sim="fixed")
  lb.MKtau <- round(sort(boot.out.MKtau$t)[(1-ci)*nsim], digits = 7)
  ub.MKtau <- round(sort(boot.out.MKtau$t)[ci*nsim], digits = 7)
  lb.MKZ <- round(sort(boot.out.Zval$t)[(1-ci)*nsim], digits = 7)
  ub.MKZ <- round(sort(boot.out.Zval$t)[ci*nsim], digits = 7)

  cat(paste("Z-Value = ", Z,
            "Sen's Slope = ", slp,
            "S = ", S,
            "Kendall's Tau = ", Tau,
            "Kendall's Tau Empirical Bootstrapped CI =", sprintf("(%s,%s)",lb.MKtau,ub.MKtau),
            "Z-value Empirical Bootstrapped CI =", sprintf("(%s,%s)",lb.MKZ,ub.MKZ),sep="\n"))

}
