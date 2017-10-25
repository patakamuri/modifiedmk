#' @title Spearman's Rank Correlation Test
#'
#' @description Spearman's Rank Correlation test by Sneyers (1990) is an useful in detecting trends.
#'
#' @param  x  - Time series data vector
#'
#' @usage spear(x)
#'
#' @return  Correlation coefficient - Spearman's Correlation coefficient value
#'
#' @return  Test Statistics - Z-Transform value to test significance (Rho*(sqrt(n-1)))
#'
#' @export
#'
#' @references Yue, S., & Wang, C. Y. (2004). The Mann-Kendall test modified by effective sample size to detect trend in serially correlated hydrological series. Water Resources Management, 18(3), 201–218. http://doi.org/10.1023/B:WARM.0000043140.61082.60
#'
#' @examples x<-c(Nile)
#' spear(x)
#'
#' @export
#'

spear<-function(x){
  x=x
# creating a sequential series whose length is equal to the input data
  xi<-1:length(x)

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
  tsrc<-rhos*sqrt(length(x))

  return(c("Correlation coefficient" = rhos, "Test Statistics" = tsrc))
}
