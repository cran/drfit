\name{drcfit}
\alias{drcfit}
\title{Fit dose-response models using the drc package}
\description{
  Fit dose-response relationships to dose-response data and calculate
  biometric results for (eco)toxicity evaluation using the drc package
}
\usage{
  drcfit(data, chooseone = TRUE, probit = TRUE, logit = FALSE,
  weibull = FALSE, linlogit = FALSE, level = 0.95,
  showED50 = FALSE, EDx = NULL)
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
  \item{probit}{
    A boolean defining if cumulative density curves of normal distributions
    are fitted against the decadic logarithm of the dose.  Default ist TRUE.
    Note that the parameter definitions used in the model are different to the
    ones used in \code{\link{drfit}}. Parameter e from \code{\link{LN.2}} is listed
    as a here, and parameter b from there is listed as b.}
  \item{logit}{
    A boolean defining if cumulative density curves of logistic distributions
    \code{\link{plogis}} are fitted to the decadic logarithm of the dose.
    Default is FALSE.
    Again the parameter definitions used in the model are different to the
    ones used in \code{\link{drfit}}. Parameter e from \code{\link{LL.2}} is listed
    as a here, and parameter b from LL.2 is listed as b.}
  \item{weibull}{
    A boolean defining if Weibull dose-response models
    (\code{\link{W1.2}} are fitted to the untransformed dose. Default is FALSE.
    Note that the results differ from the ones obtained with
    \code{\link{drfit}}, due to a different model specification.}
  \item{linlogit}{
    A boolean defining if the linear-logistic function
    \code{\link{linlogitf}} as defined by van Ewijk and Hoekstra 1993 is
    fitted to the data. Default is FALSE. Obtaining the ED50 (and EDx values
    in general) uses \code{\link{ED}} internally and does not always give a
    result.
  }
  \item{level}{
    The level for the confidence interval listed for the log ED50.}
  \item{chooseone}{
    If TRUE (default), the models are tried in the order linlogit, probit,
    logit, weibull, and the first model that produces a valid fit is used.
    If FALSE, all models that are set to TRUE and that can be fitted will be
    reported.}
  \item{EDx}{
    A vector of inhibition values x in percent for which the corresponding doses
    EDx should be reported.
  }
  \item{showED50}{
    If set to TRUE, the ED50 and its confidence interval on the original dose
    scale (not log scale) is included in the output.
  }
}
\value{
  A dataframe with the attribute \code{models} holding a list of the fitted
  dose-response models of class \code{\link{nls}}. The dataframe has at least
  one line for each substance.

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
    by \code{\link{summary.drc}}.
  }
  \item{a}{
    For the linlogit model, this is the parameter e from \code{\link{BC.4}}.
    For the probit and the logit model, this is the ED50. For the weibull
    model, this is parameter e from \code{\link{W1.2}}. Note that the Weibull
    model is fitted to the untransformed data.
  }
  \item{b}{
    For the linlogit, probit, logit and weibull models, these are the
    parameters b from \code{\link{BC.4}}, \code{\link{LN.2}},
    \code{\link{LL.2}} and \code{\link{W1.2}}, respectively.
    Note that the parameter definitions (and in the case of Weibull, the model
    used) are different to the ones used in \code{\link{drfit}}.
  }
  \item{c}{
    Only the \dQuote{linlogit} fit produces a third parameter \code{c}, which is
    the parameter f from the \code{\link{BC.4}} function.
  }
  If the parameter \code{showED50} was set to TRUE, the ED50 values and their
  confidence intervals are also included on the original dose scale.

  If one or more response leves were specified in the argument \code{EDx},
  the corresponding dose levels are given in addition, together with their
  confidence intervals as calculated by \code{\link{ED}} from the drc package.
}
\examples{
data(antifoul)
r <- drcfit(antifoul, showED50 = TRUE, EDx = c(5, 10, 20))
format(r, digits = 2)
}
\note{There is a demo for each dataset that can be accessed by
  \code{demo(dataset)}}
\seealso{
  Further examples are given in help pages to the datasets
  \code{\link{antifoul}}, \code{\link{IM1xIPC81}} and
  \code{\link{IM1xVibrio}}.
}
\author{
  Johannes Ranke \email{jranke@uni-bremen.de}
  The functionality of the drc package used under the hood in this function
  was written by Christian Ritz.
}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
