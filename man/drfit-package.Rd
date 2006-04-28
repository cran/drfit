\name{drfit-package}
\alias{drfit-package}
\docType{package}
\title{
Dose-response data evaluation
}
\description{
See \url{../DESCRIPTION}
}
\details{
There is a preliminary version of an introductory article located in \url{../doc/drfit-Rnews.pdf}.
}
\author{
Author and Maintainer: Johannes Ranke <jranke@uni-bremen.de>
}
\note{There is a demo for each dataset that can be accessed by
  \code{demo(dataset)}} 
\keyword{ package }
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\seealso{
On CRAN, there is another, more statistically sophisticated package with
similar functionality called \code{\link[drc:drc-package]{drc}}. I think 
the advantage of my package is its user-friendliness.
}
\examples{
data(antifoul)
r <- drfit(antifoul)
format(r,digits=2)
\dontrun{drplot(r,antifoul,overlay=TRUE,bw=FALSE)}
}