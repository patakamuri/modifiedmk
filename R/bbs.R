#' @title Block bootstrapped Mann-Kendall and Spearman's Rank Correlation Test
#'
#' @description Block lengths are automatically calculated using the lag of the smallest signigicant auto-correlation coefficient.  And eta value of 1 is used as a default value following the results of Khaliq et al. (2009), 2000 bootstrapped replicates are recommended following the results of Svensson et al. (2005).
#'
#' @importFrom stats acf qnorm
#'
#' @importFrom boot tsboot boot.ci
#'
#' @usage bbs(x,ci=0.95,method=c("mk","sr"),nsim=2000,eta=1)
#'
#' @param  x  - Time series data vector
#'
#' @param  ci - Confidence Interval
#'
#' @param  nsim - Number of simulations
#'
#' @param  eta - Added to the block length
#'
#' @return  Tau  - Mann-Kendall's Tau
#'
#' Confidence Interval  - Empirical bootstrapped confidence interval
#'
#' @references Box, G. E. P. and Jenkins, G. M. (1970). Time Series Analysis Forecasting and Control. Holden-Day, San Fransisco, California, 712 pp.
#'
#' @references Kendall, M. (1975). Multivariate analysis. Charles Griffin. Londres. 0-85264-234-2.
#'
#' @references Khaliq, M. N., Ouarda, T. B. M. J., Gachon, P., Sushama, L., and St-Hilaire, A. (2009). Identification of hydrological trends in the presence of serial and cross correlations: A review of selected methods and their application to annual flow regimes of Canadian rivers. Journal of Hydrology, 368: 117-130.
#'
#' @references Kundzewicz, Z. W. and Robson, A. J. (2000). Detecting Trend and Other Changes in Meteorological Data. World Climate Program – Data and Monitoring. World Meteorological Organization, Geneva  (WMO/TD-No. 1013).
#'
#' @references Kundzewicz, Z. W. and Robson, A. J. (2004). Change detection in hydrological records – a review of the methodology. Hydrological Sciences Journal, 49(1): 7-19.
#'
#' @references Lehmann, E. L. (1975). Nonparametrics: statistical methods based on ranks. Holden-Day, Inc., California, 457 pp.
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3), 245–259. <doi:10.1017/CBO9781107415324.004>
#'
#' @references Önöz, B. and Bayazit M. (2012). Block bootstrap for Mann-Kendall trend test of serially dependent data. Hydrological Processes, 26: 3552-3560.
#'
#' @references Sneyers, R. (1990). On the statistical analysis of series of observations. World Meteorological Organization, Technical Note no. 143, WMO no. 415, 192 pp.
#'
#' @references Svensson, C., Kundzewicz, Z. W., and Maurer, T. (2005). Trend detection in river flow series: 2. Floods and low-flow index series. Hydrological Sciences Journal, 50(5): 811-823.
#'
#' @details The block bootstrap is used along with the non-parametric Mann-Kendall trend test and Spearman's Rank Correlation Test.  A test statistic falling in the tails of the simulated empirical distribution, the results is likely significant. The default technique is the Mann-Kendall trend test.
#'
#' @examples x<-c(Nile)
#' bbs(x)
#'
#' @export
#'
bbs <- function(x,ci=0.95,method=c("sr"),nsim=40,eta=1) {
  # Initialize the test Parameters

  # Time-Series Vector
  x = x
  # Confidance Interval
  ci = ci
  #Number of Simulations
  nsim=nsim
  #Value of eta
  eta=eta
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

  # To test whether the data values are finite numbers and attempting to eliminate non-finite numbers

  if (any(is.finite(x) == FALSE)) {
    x[-c(which(is.finite(x) == FALSE))] -> x
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }

  #bounds of the confidence intervals of the acf function

  bd <- qnorm((1 + ci)/2)/sqrt(n)

  # Calculating auto-correlation function of the ranks of observations (ro)

  ro <- acf(x, lag.max=round(n/4), plot=T)$acf[-1]

  #Initialize vector of significant lags of auto-correlation

  sig.v <- rep(0, round(n/4))

  #Identify the lag of the smallest significant auto-correlation coefficient

  for (i in 1:round(n/4)) {
    if (-bd > ro[i] | bd < ro[i]) {
      sig.v[i]<-ro[i]
    }
  }

  if (length(sig.v) == 1) {
    min.sig<-0
  } else {
    min.sig <- max(which(sig.v %in% min(sig.v[sig.v > 0])))
  }

  #Block length

  bl.len <- min.sig + eta

  #Block bootstrap using Mann Kendall

  if (method == "mk") {
    MK.orig <- mkttest(x)
    z<-round(MK.orig["Z-Value"], digits = 7)
    slp<-round(MK.orig["Sen's slope"], digits = 7)
    S<-MK.orig["S"]
    MKtau <- function(x) mkttest(x)[["Tau"]]
    boot.out <- tsboot(x, MKtau, R=nsim, l=bl.len, sim="fixed")
    Tau <- round(boot.out$t0, digits = 7)
    bbs.ci <- boot.ci(boot.out, conf = ci, type="basic")$basic[4:5]
    lb <- round(bbs.ci[1], digits = 7)
    ub <- round(bbs.ci[2], digits = 7)

    cat(paste("Z-Value = ", z,
              "Sen's slope = ", slp,
              "S = ", S,
              "Kendall's Tau = ", Tau,
              "Kendall's Tau Empirical Bootstrapped CI =", sprintf("(%s,%s)",lb,ub),sep="\n"),sep="\n")

  } else if (method == "sr") {
    SR.orig<-round(spear(x)["Correlation coefficient"], digits = 7)
    SRfunc <- function(x) spear(x)[["Test Statistics"]]
    boot.out <- tsboot(x, SRfunc, R=nsim, l=bl.len, sim="fixed")
    SRtsrc <- round(boot.out$t0, digits = 7)
    bbs.ci <- boot.ci(boot.out, conf = ci, type="basic")$basic[4:5]
    lb <- round(bbs.ci[1], digits = 7)
    ub <- round(bbs.ci[2], digits = 7)

    cat(paste("Spearman's Correlation Coefficient = ", SR.orig,
              "Test Statistic = ", SRtsrc ,
              "Test Statistic Empirical Bootstrapped CI =", sprintf("(%s,%s)",lb,ub),sep="\n"),sep="\n")
  }
}

