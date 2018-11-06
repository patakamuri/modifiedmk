#' @title Hamed (2009) Bias Corrected Prewhitening.
#'
#' @description Hamed (2009) proposed a prewhitening technique in which the slope and lag-1 serial correltaion coefficient are simultaneously estimated.  The lag-1 serial correltaion coefficient is then corrected for bias before prewhitening.
#'
#' @importFrom utils head tail
#'
#' @importFrom stats pnorm qnorm
#'
#' @usage bcpw(x)
#'
#' @param  x  - Time series data vector
#'
#' @return  Z-Value - Mann-Kendall Z-statistic after bias corrected prewhitening
#'
#' @return  Prewhitened Sen's Slope - Sen's slope of the prewhitened data
#'
#' @return  Sen's Slope - Sen's slope for the original data series 'x'
#'
#' @return  P-value - p-value after prewhitening
#'
#' @return  S - Mann-Kendall 'S' statistic
#'
#' @return  Var(s) - Variance of 'S'
#'
#' @return  Tau - Mann-Kendall's Tau
#'
#' @references Hamed, K. H. (2009). Enhancing the effectiveness of prewhitening in trend analysis of hydrologic data. Journal of Hydrology, 368: 143-155.
#'
#' @references Kendall, M. (1975). Multivariate analysis. Charles Griffin. Londres. 0-85264-234-2.
#'
#' @references Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3), 245-259. <doi:10.1017/CBO9781107415324.004>
#'
#' @references van Giersbergen, N. P. A. (2005). On the effect of deterministic terms on the bias in stable AR models. Economic Letters, 89: 75-82.
#'
#' @details Employs ordinary least squares (OLS) to simultaneously estimate the lag-1 serial correlation coefficient and slope of trend.  The lag-1 serial correlation coefficient is then bias corrected.
#'
#' @examples x<-c(Nile)
#' bcpw(x)
#'
#' @export
#'
bcpw <- function(x) {
  # Initialize the test Parameters

  # Time-Series Vector
  x = x
  # Modified Z-Statistic after Pre-Whitening
  z = NULL
  # Modified P-value after Pre-Whitening
  pval = NULL
  # Initialize Mann-Kendall 'S'- Statistic
  S = 0
  # Initialize Mann-Kendall var.S
  var.S = NULL
  # Initialize Mann-Kendall Tau
  Tau = NULL

  # To test whether the data is in vector format

  if (is.vector(x) == FALSE) {
    stop("Input data must be a vector")
  }

  nx<-length(x)

  # To test whether the data values are finite numbers and attempting to eliminate non-finite numbers
  if (any(is.finite(x) == FALSE)) {
    x[-c(which(is.finite(x) == FALSE))] -> x
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }

  #Calculate the lag 1 autocorrelation coefficient and the intercept
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

  # Calculating pre-whitened Series

  a=1:(nx-1)
  b=2:nx
  xn<-(x[b]-(x[a]*ACFlag1BC))

  PWn<-length(xn)

  # Calculating Mann-Kendall 'S'- Statistic

  for (i in 1:(PWn-1)) {
    for (j in (i+1):PWn) {
      S = S + sign(xn[j]-xn[i])
    }
  }

  # Calculating Mann-Kendall Variance (Var(s))

  var.S = PWn*(PWn-1)*(2*PWn+5)*(1/18)
  if(length(unique(xn)) < PWn) {
    unique(xn) -> aux
    for (i in 1:length(aux)) {
      length(which(xn == aux[i])) -> tie
      if (tie > 1) {
        var.S = var.S - tie*(tie-1)*(2*tie+5)*(1/18)
      }
    }
  }

  # Calculating Z-Statistic values

  if (S == 0) {
    z = 0
  }
  if (S > 0) {
    z = (S-1)/sqrt(var.S)
  } else {
    z = (S+1)/sqrt(var.S)
  }

  # Calculating P-Value before and after prewhitening

  pval = 2*pnorm(-abs(z))


  # Calculating kendall's Tau

  Tau = S/(.5*PWn*(PWn-1))


  # Calculating Sen's slope for original series 'x'

  rep(NA, nx * (nx - 1)/2) -> V
  k = 0
  for (i in 1:(nx-1)) {
    for (j in (i+1):nx) {
      k = k+1
      V[k] = (x[j]-x[i])/(j-i)
    }
  }
  median(V,na.rm=TRUE)->slp

  # Calculating Sen's slope for PW series

  rep(NA, PWn * (PWn - 1)/2) -> W
  m = 0
  for (i in 1:(PWn-1)) {
    for (j in (i+1):PWn) {
      m = m+1
      W[m] = (xn[j]-xn[i])/(j-i)
    }
  }
  median(W,na.rm=TRUE)->slp1

  cat(paste("Z-Value = ", z,
            "Prewhitened Sen's Slope = ", slp1,
            "Sen's Slope = ", slp,
            "P-value = ", pval,
            "S = ", S,
            "Var(S) = ", var.S,
            "Tau = ", Tau))
}

