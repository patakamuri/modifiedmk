#' @title Bootstrapped Mann-Kendall Trend Test with Optional Bias Corrected Prewhitening
#'
#' @description The empirical distribution of the Mann-Kendall test statistic is calculated by bootstrapped resampling.  The Hamed (2009) bias correction prewhitening technique can optionally be applied as the default for prewhitening before the bootstrapped Mann-Kendall test is applied (Lacombe et al., 2012).
#'
#' @importFrom boot tsboot
#'
#' @usage pbmk(x, nsim=1000, pw="Hamed")
#'
#' @param  x  - Time series data vector
#'
#' @param  nsim - Number of bootstrapped simulations
#'
#' @param pw -  Optional bias corrected prewhitening suggested by Hamed (2009)
#'
#' @return  Z Value - Mann-Kendall Z statistic from original data
#' 
#' @return  Sen's Slope - Sen's slope from the original data
#' 
#' @return  S - Mann-Kendall S statistic
#' 
#' @return  Kendall's Tau - Mann-Kendall's Tau
#' 
#' @return  BCP Z Value - Bias corrected prewhitened Z value
#' 
#' @return  BCP Sen's Slope - Bias corrected prewhitened Sen's slope
#' 
#' @return  BCP S - Bias corrected prewhitened S
#' 
#' @return  BCP Kendall's Tau - Bias corrected prewhitened Kendall's Tau
#' 
#' @return  Bootstrapped P-Value - Mann-Kendall bootstrapped p-value
#'
#' @references Hamed, K. H. (2009). Enhancing the effectiveness of prewhitening in trend analysis of hydrologic data. Journal of Hydrology, 368: 143-155.
#'
#' @references Kendall, M. (1975). Rank Correlation Methods. Griffin, London, 202 pp.
#'
#' @references Kundzewicz, Z. W. and Robson, A. J. (2004). Change detection in hydrological records - a review of the methodology. Hydrological Sciences Journal, 49(1): 7-19.
#'
#' @references Lancombe, G., McCartney, M., and Forkuor, G. (2012). Drying climate in Ghana over the period 1960-2005: evidence from the resampling-based Mann-Kendall test at local and regional levels. Hydrological Sciences Journal, 57(8): 1594-1609.
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3): 245-259.
#'
#' @references van Giersbergen, N. P. A. (2005). On the effect of deterministic terms on the bias in stable AR models. Economic Letters, 89: 75-82.
#'
#' @references Yue, S. and Pilon, P. (2004). A comparison of the power of the t test, Mann-Kendall and bootstrap tests for trend detection, Hydrological Sciences Journal, 49(1): 21-37.
#'
#' @details Bootstrapped samples are calculated by resampling one value at a time from the time series with replacement.  The p-value (\eqn{p_s}) of the resampled data is estimated by (Yue and Pilon, 2004): \deqn{p_s = m_s/M} The Mann-Kendall test statistics (S) is calculated for each resampled dataset.  The resultant vector of resampled S statistics is then sorted in ascending ordering, where \eqn{p_s} is the rank corresponding the largest bootstrapped value of S being less than the test statistic value calculated from the actual data.  M is the total number of bootstrapped resamples.  The default value of M is 1000, however, Yue and Pilon (2004) suggest values between 1000 and 2000. If the user does not choose to apply prewhitening, this argument 'pw' can be set to NULL.
#' 
#'
#' @examples x<-c(Nile[1:10])
#' pbmk(x)
#'
#' @export
#'
pbmk <- function(x, nsim=1000, pw="Hamed") {
  # Initialize the test parameters

  # Time series vector
  x = x
  #Number of simulations
  nsim=nsim
  # Mann-Kendall Tau
  Tau = NULL
  # Modified Z statistic after prewhitening
  Z = NULL
  # Modified p-value after prewhitening
  pval = NULL
  # Initialize Mann-Kendall S statistic - prewhitened
  S = NULL
  # Initialize Mann-Kendall S statistic
  S.orig = NULL
  # Sen's slope estimate
  slp = NULL

  # To test whether the data is in vector format

  if (is.vector(x) == FALSE) {
    stop("Input data must be a vector")
  }

  nx<-length(x)

  #Specify minimum input vector length
  if (nx < 3) {
    stop("Input vector must contain at least three values")
  }

  # To test whether the data values are finite numbers and attempting to eliminate non-finite numbers
  if (any(is.finite(x) == FALSE)) {
    x[-c(which(is.finite(x) == FALSE))] -> x
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }
  
  nx<-length(x)

  if (is.null(pw) == FALSE) {
    #Calculate the lag-1 autocorrelation coefficient and the intercept
    zx<-cbind(head(x,n=nx-1),matrix(data=1, nrow=(nx-1),ncol=1),tail(seq(1:nx),n=(nx-1)))
    y<-tail(x,n=nx-1)
    zTrans<-t(zx)
    zTransz<-zTrans%*%zx
    zTranszInv<-solve(zTransz)
    zTranszInvzTrans<-zTranszInv%*%zTrans
    params<-zTranszInvzTrans%*%y
    ACFlag1<-params[1]

    #Correct for bias in the lag-1 acf using eq. 24 of Hamed (2009)
    ACFlag1BC<-((nx*ACFlag1)+2)/(nx-4)

    # Calculating prewhitened series

    a=1:(nx-1)
    b=2:nx
    xn<-(x[b]-(x[a]*ACFlag1BC))

    #Bootstrapped using Mann-Kendall
    MK.orig <- mkttest(x)
    Z <- round(MK.orig["Z-Value"], digits = 7)
    slp <- round(MK.orig["Sen's slope"], digits = 7)
    Tau <- round(MK.orig["Tau"], digits = 7)
    S.orig <- MK.orig["S"]
    MKpw <- mkttest(xn)
    Zpw <- round(MKpw["Z-Value"], digits = 7)
    slpPW <- round(MKpw["Sen's slope"], digits = 7)
    TauPW <- round(MKpw["Tau"], digits = 7)
    Spw <- MKpw["S"]
    MKS <- function(xn) mkttest(xn)[["S"]]
    boot.out.MKS <- tsboot(xn, MKS, R=nsim, l=1, sim="fixed")
    loc <- suppressWarnings(max(which(sort(boot.out.MKS$t) < Spw)))
    if (loc == -Inf) {
      loc <- 1
    }
    pval <- loc/nsim

    cat(paste("Z Value = ", Z,
              "Sen's Slope = ", slp,
              "S = ", S.orig,
              "Kendall's Tau = ", Tau,
              "BCP Z Value = ", Zpw,
              "BCP Sen's Slope = ", slpPW,
              "BCP S = ", Spw,
              "BCP Kendall's Tau = ", TauPW,
              "Bootstrapped P-Value =", pval ,sep="\n"))
  } else {
    #Bootstrapped using Mann-Kendall
    MK.orig <- mkttest(x)
    Z <- round(MK.orig["Z-Value"], digits = 7)
    slp <- round(MK.orig["Sen's slope"], digits = 7)
    Tau <- round(MK.orig["Tau"], digits = 7)
    S.orig <- MK.orig["S"]
    MKS1 <- function(x) mkttest(x)[["S"]]
    boot.out.MKS1 <- tsboot(x, MKS1, R=nsim, l=1, sim="fixed")
    loc <- suppressWarnings(max(which(sort(boot.out.MKS1$t) < S.orig)))

    if (loc == -Inf) {
      loc <- 1
    }
    pval <- loc/nsim

    cat(paste("Z Value = ", Z,
              "Sen's Slope = ", slp,
              "S = ", S.orig,
              "Kendall's Tau = ", Tau,
              "Bootstrapped P-Value =", pval ,sep="\n"))
  }
}



