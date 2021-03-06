% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmkh.R
\name{mmkh}
\alias{mmkh}
\title{Modified Mann-Kendall Test For Serially Correlated Data Using the Hamed and Rao (1998) Variance Correction Approach}
\usage{
mmkh(x,ci=0.95)
}
\arguments{
\item{x}{- Time series data vector}

\item{ci}{- Confidence interval}
}
\value{
Corrected Zc  - Z statistic after variance Correction

new P.value  - P-value after variance correction

N/N*  - Effective sample size

Original Z  - Original Mann-Kendall Z statistic

Old P-value  - Original Mann-Kendall p-value

Tau  - Mann-Kendall's Tau

Sen's Slope  - Sen's slope

old.variance - Old variance before variance Correction

new.variance - Variance after correction
}
\description{
Time series data is often influenced by previous observations. When data is not random and influenced by autocorrelation, modified Mann-Kendall tests may be used for trend detction studies. Hamed and Rao (1998) have proposed a variance correction approach to address the issue of serial correlation in trend analysis. Data are initially detrended and the effective sample size is calulated using the ranks of significant serial correlation coefficients which are then used to correct the inflated (or deflated) variance of the test statistic.
}
\details{
A detrended time series is constructed using Sen's slope and the lag-1 autocorreltation coefficient of the ranks of the data. The variance correction approach proposed by Hamed and Rao (1998) uses only significant lags of autocorrelation coefficients.
}
\examples{
x<-c(Nile)
mmkh(x)

}
\references{
Hamed, K. H. and Rao, A. R. (1998). A modified Mann-Kendall trend test for autocorrelated data. Journal of Hydrology, 204(1–4): 182–196. <doi:10.1016/S0022-1694(97)00125-X>

Kendall, M. (1975). Rank Correlation Methods. Griffin, London, 202 pp.

Mann, H. B. (1945). Nonparametric Tests Against Trend. Econometrica, 13(3): 245-259.

Rao, A. R., Hamed, K. H., & Chen, H.-L. (2003). Nonstationarities in hydrologic and environmental time series. Ringgold Inc., Portland, Oregon, 362 pp. <doi:10.1007/978-94-010-0117-5>

Salas, J.D. (1980). Applied modeling of hydrologic times series. Water Resources Publication, 484 pp.

Sen, P. K. (1968). Estimates of the Regression Coefficient Based on Kendall’s Tau. Journal of the American statistical Association, 63(324): 1379. <doi:10.2307/2285891>
}
