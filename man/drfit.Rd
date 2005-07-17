\name{drfit}
\alias{drfit}
\title{Fit dose-response models}
\description{
	Fit dose-response relationships to dose-response data and calculate
  biometric results for (eco)toxicity evaluation
}
\usage{
  drfit(data, startlogEC50 = NA, chooseone = TRUE, probit = TRUE, logit = FALSE,
    weibull = FALSE, linlogit = FALSE, linlogitWrong = NA, allWrong = NA, 
    s0 = 0.5, b0 = 2, f0 = 0)
}
\arguments{
  \item{data}{
    A data frame containing dose-response data. The data frame has to
    contain at least a factor called "substance", a vector called "unit"
    containing the unit used for the dose, a column "response" with the
    response values of the test system normalized between 0 and 1 and a column
    "dose" with the numeric dose values. Such a data frame can be easily obtained
    if a compliant RODBC data source is available for use in conjunction with 
    the function \code{\link{drdata}}.  
    
    If there is a column called "ok" and it is set to "no fit" in a specific
    line, then the corresponding data point will be excluded from the fitting
    procedure, although it will be plotted.
    }
  \item{startlogEC50}{
    Especially for the linlogit model, a suitable log10 of the EC50 has to be given 
    by the user, since it is not correctly estimated for data showing hormesis with
    the default estimation method.}
  \item{probit}{
    A boolean defining if cumulative density curves of normal distributions
    \code{\link{pnorm}} are fitted against the decadic logarithm of the dose.
    Default ist TRUE.} 
  \item{logit}{
    A boolean defining if cumulative density curves of logistic distributions
    \code{\link{plogis}} are fitted to the decadic logarithm of the dose.
    Default is FALSE.} 
  \item{weibull}{
    A boolean defining if the cumulative density curves of weibull distributions
    (\code{\link{pweibull}} with additionall location parameter and scale=1)
    are fitted to the decadic logarithm of the dose. Default is FALSE.}
  \item{linlogit}{
    A boolean defining if the linear-logistic function
    \code{\link{linlogitf}} as defined by van Ewijk and Hoekstra 1993 is
    fitted to the data. Default is FALSE.}
  \item{linlogitWrong}{
    An optional vector containing the names of the substances for which the
    linlogit function produces a wrong fit.}
  \item{allWrong}{
    An optional vector containing the names of the substances for which all
    functions produce a wrong fit.}
  \item{chooseone}{
    If TRUE (default), the models are tried in the order linlogit, weibull,
    logit, probit, and the first model that produces a valid fit is used.
    If FALSE, all models that are set to TRUE and that can be fitted will be
    reported.}
  \item{s0}{
    If the weibull model is fitted, s0 gives the possibility to adjust the
    starting value for the shape parameter of \code{\link{pweibull}}.}
  \item{b0,f0}{
    If the linearlogistic model is fitted, b0 and f0 give the possibility to
    adjust the starting values for the parameters b and f.}
}
\value{
  \item{results}{
      A data frame containing at least one line for each substance. If the data
      did not show a mean response < 0.5 at the highest dose level, the
      modeltype is set to "none".
      Every successful fit is reported in one line. Parameters of the fitted
      curves are only reported if the fitted EC50 is not higher than the
      highest dose. 
      \code{n} is the number of dose-response curves in the raw data (repetitions
      in each point), \code{lld} is the decadic logarithm of the lowest dose and
      \code{lhd} is the decadic logarithm of the highest dose.
      For the "linlogit", "logit" and "probit" models, the parameter
      \code{a} that is reported coincides with the logEC50, i.e the logEC50 is 
      one of the model parameters that is being fitted, and therefore
      a standard deviation \code{std} is reported for the logEC50. In the 
      case of the "weibull" model, \code{a} is a location parameter.
      Parameter \code{b} in the case of the "linlogit" fit is the variable
      b from the \code{\link{linlogitf}} function. In the case of "probit" fit
      it is the standard deviation of the fitted normal distribution, in the case
      of the "logit" fit it is the \code{scale} parameter in the \code{\link{plogis}}
      function, and in the "weibull" fit it is the \code{shape} parameter of the
      fitted \code{\link{pweibull}} function. Only the "linlogit" fit produces a 
      third parameter \code{c} which is the variable f from the
      \code{\link{linlogitf}} function.}
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
