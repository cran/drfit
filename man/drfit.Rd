\name{drfit}
\alias{drfit}
\title{Fit dose-response models}
\description{
  Fit dose-response relationships to dose-response data and calculate
  biometric results for (eco)toxicity evaluation
}
\usage{
  drfit(data, startlogED50 = NA, chooseone = TRUE, probit = TRUE, logit = FALSE,
  weibull = FALSE, linlogit = FALSE, level = 0.95, linlogitWrong = NA,
  allWrong = NA, ps0 = 1, ls0 = 0.5, ws0 = 0.5, b0 = 2, f0 = 0,
  showED50 = FALSE,
  EDx = NULL, EDx.tolerance = 1e-4)
}
\arguments{
  \item{data}{
    A data frame containing dose-response data. The data frame has to contain
    at least a factor called \dQuote{substance}, a numeric vector \dQuote{dose}
    with the dose values, a vector called \dQuote{unit} containing the unit
    used for the dose and a numeric vector \dQuote{response} with the response
    values of the test system normalized between 0 and 1. Such a data frame can
    be easily obtained if a compliant RODBC data source is available for use in
    conjunction with the function \code{\link{drdata}}.

    If there is a column called \dQuote{ok} and it is set to \dQuote{no fit} in
    a specific line, then the corresponding data point will be excluded from
    the fitting procedure, although it will be plotted.}
  \item{startlogED50}{
    Especially for the linlogit model, a suitable log10 of the ED50 has to be
    given by the user, since it is not correctly estimated for data showing
    hormesis with the default estimation method.}
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
    are fitted to the decadic logarithm of the dose. Default is FALSE.
    Note that the weibull distribution is fitted here to the log transformed doses
    which appears to be an uncommon approach.}
  \item{linlogit}{
    A boolean defining if the linear-logistic function
    \code{\link{linlogitf}} as defined by van Ewijk and Hoekstra 1993 is
    fitted to the data. Default is FALSE.}
  \item{level}{
    The level for the confidence interval listed for the log ED50.}
  \item{linlogitWrong}{
    An optional vector containing the names of the substances for which the
    linlogit function produces a wrong fit.}
  \item{allWrong}{
    An optional vector containing the names of the substances for which all
    functions produce a wrong fit.}
  \item{chooseone}{
    If TRUE (default), the models are tried in the order linlogit, probit,
    logit, weibull, and the first model that produces a valid fit is used.
    If FALSE, all models that are set to TRUE and that can be fitted will be
    reported.}
  \item{ps0}{
    If the probit model is fitted, \code{ps0} gives the possibility to adjust
    the starting value for the scale parameter of \code{\link{pnorm}}.}
  \item{ls0}{
    If the logit model is fitted, \code{ls0} gives the possibility to adjust
    the starting value for the scale parameter of \code{\link{plogis}}.}
  \item{ws0}{
    If the weibull model is fitted, \code{ws0} gives the possibility to adjust
    the starting value for the shape parameter of \code{\link{pweibull}}.}
  \item{b0,f0}{
    If the linearlogistic model is fitted, \code{b0} and \code{f0} give the
    possibility to adjust the starting values for the parameters b and f.}
  \item{showED50}{
    If set to TRUE, the ED50 and its confidence interval on the original dose
    scale (not log scale) is included in the output.
  }
  \item{EDx}{
    A vector of inhibition values x in percent for which the corresponding doses
    EDx should be reported.
  }
  \item{EDx.tolerance}{
    Tolerance of the effect level, expressed on the response scale from 0 to 1.
  }
}
\value{
  A dataframe with the attribute \code{models} holding a list of the fitted
  dose-response models of class \code{\link{nls}}. The dataframe has at least
  one line for each substance.

  For the \dQuote{linlogit}, \dQuote{logit} and \dQuote{probit} models, the
  parameter \code{a} that is reported coincides with the logED50, i.e the
  logED50 is one of the model parameters that is being fitted. Therefore,
  a confidence interval for the confidence level \code{level} is calculated
  using the \code{\link[MASS:confint]{confint.nls}} function and listed.

  The following variables are in the dataframe:
  \item{Substance}{
    The name of the substance
  }
  \item{ndl}{
    The number of dose levels in the raw data
  }
  \item{n}{
    The total number of data points in the raw data used for the fit
  }
  \item{lld}{
    The decadic logarithm of the lowest dose
  }
  \item{lhd}{
    The total number of data points in the raw data used for the fit
  }
  \item{mtype}{
    If the data did not show a mean response < 0.5 at the highest dose level,
    the modeltype is set to \dQuote{inactive}. If the mean response at the
    lowest dose is smaller than 0.5, the modeltype is set to \dQuote{active}.
    In both cases, no fitting procedure is carried out. If the fitted ED50
    is higher than the highest dose, \dQuote{no fit} is given here.
  }
  \item{logED50}{
    The decadic logarithm of the ED50
  }
  \item{low \%}{
    The lower bound of the confidence interval of log ED50.
    The name of the column depends on the requested confidence \code{level}.
  }
  \item{high \%}{
    The higher bound of the confidence interval of log ED50.
    The name of the column depends on the requested confidence \code{level}.
  }
  \item{unit}{
    The unit used for the dose levels in the dose-response data
  }
  \item{sigma}{
    The square root of the estimated variance of the random error as returned
    by \code{\link{summary.nls}}.
  }
  \item{a}{
    For the \dQuote{linlogit}, \dQuote{logit} and \dQuote{probit} models, the
    parameter \code{a} coincides with the logED50.  In the case of the
    \dQuote{weibull} model, \code{a} is a location parameter.
  }
  \item{b}{
    Parameter \code{b} in the case of the \dQuote{linlogit} fit is the variable
    b from the \code{\link{linlogitf}} function. In the case of \dQuote{probit}
    fit it is the standard deviation of the fitted normal distribution, in the
    case of the \dQuote{logit} fit it is the \code{scale} parameter in the
    \code{\link{plogis}} function, and in the \dQuote{weibull} fit it is the
    \code{shape} parameter of the fitted \code{\link{pweibull}} function.
  }
  \item{c}{
    Only the \dQuote{linlogit} fit produces a third parameter \code{c} which is
    the variable f from the \code{\link{linlogitf}} function.
  }
  If the parameter \code{showED50} was set to TRUE, the ED50 values and their
  confidence intervals are also included on the original dose scale.

  If one or more response leves were specified in the argument \code{EDx},
  the corresponding dose levels are given in addition.
}
\examples{
data(antifoul)
r <- drfit(antifoul, showED50 = TRUE, EDx = c(5, 10, 20))
format(r, digits = 2)
}
\note{There is a demo for each dataset that can be accessed by
  \code{demo(dataset)}}
\seealso{
  Further examples are given in help pages to the datasets
  \code{\link{antifoul}}, \code{\link{IM1xIPC81}} and
  \code{\link{IM1xVibrio}}.
  Since version 0.6.1 of this package, there is a drop-in replacement function
  \code{\link{drcfit}} which internally uses the drc package and also gives
  confidence intervals for EDx values via this package.
}
\author{
  Johannes Ranke
  \email{jranke@uni-bremen.de}
  \url{http://www.uft.uni-bremen.de/chemie/ranke}
}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
