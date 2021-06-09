#' @title Spearman's Rank Correlation Test
#'
#' @description Spearman's Rank Correlation test by Lehmann (1975) and Sneyers (1990) is useful in detecting trends.
#'
#' @param  x  - Time series data vector
#'
#' @usage spear(x)
#'
#' @return  Correlation coefficient - Spearman's Correlation coefficient value
#'
#' @return Z-Tranformed Test Statistic value - Z-transform value to test significance \eqn{\rho(\sqrt{n-1})}
#'
#' @export
#'
#' @references Lehmann, E. L. (1975). Nonparametrics: statistical methods based on ranks. Holden-Day, Inc., California, 457 pp.
#'
#' @references Sneyers, R. (1990). On the statistical analysis of series of observations. World Meteorological Organization, Technical Note no. 143, WMO no. 415, 192 pp.
#'
#' @details Spearman's Rank Correlation test by Lehmann (1975) and Sneyers (1990) is implemeted in this function.
#'
#' @examples x<-c(Nile)
#' spear(x)
#'
#' @export
#'

spear<-function(x) {
  # Initialize the test parameters
  options(scipen = 999)
  # Time series vector
  x=x
  #Correlation coefficient
  rhos=NULL
  #Z-tranformed Test Statistic value
  tsrc=NULL

  # creating a sequential series whose length is equal to the input data
  xi<-1:length(x)

  # To test whether the data is in vector format

  if (is.vector(x) == FALSE) {
    stop("Input data must be a vector")
  }

  # To test whether the data values are finite numbers and attempting to eliminate non-finite numbers

  if (any(is.finite(x) == FALSE)) {
    x[-c(which(is.finite(x) == FALSE))] -> x
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }

  n<-length(x)

  #Specify minimum input vector length
  if (n < 3) {
    stop("Input vector must contain at least three values")
  }

  #Calculating ranks of the data
  yi<-as.integer(rank(x))

  # calculating Sx,Sy and Sxy

  for (i in 1: length(xi)){
    sx=0
    sx=sx+((xi-mean(xi))^2)
  }


  for (i in 1: length(yi)){
    sy=0
    sy=sy+((yi-mean(yi))^2)
  }

  for (i in 1: length(yi)){
    sxy=0
    sxy=sxy+((xi-mean(xi))*(yi-mean(yi)))
  }

 # Calculating Rho(s)

  rhos<- sum(sxy)/sqrt(sum(sx)*sum(sy))

 # Calculating Rho(s)*sqrt(n-1)
  tsrc<-rhos*sqrt(length(x)-1)

  return(c("Correlation coefficient" = rhos,
           "Z-tranformed Test Statistic value" = tsrc))

  }
