\name{drfit}
\alias{drfit}
\title{Fit dose-response models}
\description{
	Fit dose-response relationships to dose-response data and calculate
  biometric results for (eco)toxicity evaluation
}
\usage{
  drfit(data, startlogEC50 = NA, chooseone = TRUE, lognorm = TRUE, logis = FALSE, 
    linearlogis = FALSE, linearlogisWrong = NA, b0 = 2, f0 = 0)
}
\arguments{
  \item{data}{
    A data frame as returned from \code{\link{drdata}}.  The data frame has to
    contain at least a factor called "substance", a vector called "unit"
    containing the unit used for the dose, a column "response" with the
    response values of the test system normalized between 0 and 1 and a column
    "dose" with the numeric dose values. For later use of the
    \code{\link{drplot}} function, a factor called "dosefactor" also has to be
    present, containing the dose as a factor.
    }
  \item{startlogEC50}{
    Especially for the linearlogis model, a suitable log10 of the EC50 has to be given 
    by the user, since it is not correctly estimated for data showing hormesis with
    the default estimation method.}
  \item{lognorm}{
    A boolean defining if cumulative density curves of normal distributions
    are fitted to the data. Default ist TRUE.} 
  \item{logis}{
    A boolean defining if cumulative densitiy curves of logistic distributions
    are fitted to the data. Default is FALSE.} 
  \item{linearlogis}{
    A boolean defining if the linear-logistic function
    \code{\link{linearlogisf}} as defined by van Ewijk and Hoekstra 1993 is
    fitted to the data. Default is FALSE.}
  \item{linearlogisWrong}{
    An optional vector containing the names of the substances for which the
    linearlogis function produces a wrong fit.}
  \item{chooseone}{
    If TRUE (default), the models are tried in the order linearlogis, logis and lognorm,
    and the first model that produces a valid fit is used. Usually this will be the one 
    with the lowest residual standard deviation. If FALSE, all models that are set to TRUE
    and that can be fitted will be reported.}
  \item{b0,f0}{
    If the linearlogistic model is fitted, b0 and f0 give the possibility to
    adapt the starting values for the parameters b and f.}
}
\value{
  \item{results}{
      A data frame containing at least one line for each substance. If the data did not
      show a mean response < 0.5 at the highest dose level, the modeltype is set to "none".
      Every successful fit is reported in one line. Parameters of the fitted curves are only
      reported if the fitted EC50 is not higher than the highest dose.}
    
} 
\examples{
data(antifoul)
r <- drfit(antifoul)
format(r,digits=2)
}
\author{
  Johannes Ranke 
  \email{jranke@uni-bremen.de} 
  \url{http://www.uft.uni-bremen.de/chemie/ranke}
}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
