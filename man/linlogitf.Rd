\name{linlogitf}
\alias{linlogitf}
\title{Linear-logistic function}
\description{
  Helper function describing a special type of dose-response curves, showing a stimulus
  at subtoxic doses.
}
\usage{
  linlogitf(x,k,f,mu,b)
}
\arguments{
  \item{x}{
    In this context, the x variable is the dose.}
  \item{k}{
    In the drfit functions, k is set to 1.}
  \item{f}{
    One of the parameters describing the curve shape.}
  \item{mu}{
    The parameter describing the location of the curve (log ED50).}
  \item{b}{
    One of the parameters describing the curve shape.}
}
\value{
  The response at dose x.
}
\references{
  van Ewijk, P. H. and Hoekstra, J. A. (1993) \emph{Ecotox Environ Safety}
  \bold{25} 25-32}
\author{
  Johannes Ranke
  \email{jranke@uni-bremen.de}
  \url{http://www.uft.uni-bremen.de/chemie/ranke}
}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
