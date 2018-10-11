\name{drplot}
\alias{drplot}
\title{Plot dose-response models}
\description{
  Produce graphics of dose-response data and dose-response relationships
  either combined or separately, for one or more substances.
}
\usage{
  drplot(drresults, data, dtype, alpha, ctype, path,
     fileprefix, overlay, xlim, ylim, xlab, ylab, axes, frame.plot, postscript,
     pdf, png, bw, pointsize, colors, ltys, pchs, devoff, lpos)
}
\arguments{
  \item{drresults}{
    A data frame as returned from \code{\link{drfit}}.
    }
  \item{data}{
    A data frame as returned from \code{\link{drdata}}. The data frame has to
    contain at least a factor called "substance", a vector called "unit"
    containing the unit used for the dose, a column "response" with the
    response values of the test system normalized between 0 and 1, a column
    "dose" with the numeric dose values and a factor called "dosefactor"
    containing the dose as a factor.
    }
  \item{dtype}{
    A string describing if the raw data should be plotted ("raw"), or an error
    bar should be constructed from the standard deviations of the responses at
    each dose level ("std", default value) or from the confidence intervals
    ("conf"). If you don't want to see the data, set it to "none".
    }
  \item{alpha}{
    The confidence level, defaulting to 0.95, only used if dtype "conf" has been
    chosen.
    }
  \item{ctype}{
    This argument decides if horizontal lines are drawn to show the scatter of
    the control values (dose = 0), if there are more than three of them.
    Defaults to "none", further allowed values are "std" and "conf" for
    displaying the standard deviation of the controls or the confidence
    interval for the mean of the controls.
    }
  \item{path}{
    The path where graphic files should be put if any are produced. Defaults
    to "./" i.e. the current working directory of R.
    }
  \item{fileprefix}{
    A string which will form the beginning of each filename, if graphic files
    are created. Defaults to "drplot".
    }
  \item{overlay}{
    If TRUE, all output will be put into one graph, otherwise a separate graph
    will be created for each substance.
    }
  \item{xlim}{
    The plot limits (min,max) on the dose axis.
  }
  \item{ylim}{
    The plot limits (min,max) on the response axis.
  }
  \item{xlab}{
    The axis title for the x axis. Defaults to "Concentration in" \code{unit}.
  }
  \item{ylab}{
    The axis title for the y axis. Defaults to "Normalized response".
  }
  \item{axes}{
    Specifies if axes should be drawn. Default is TRUE, as in
    \code{\link{plot.default}}
  }
  \item{frame.plot}{
    Specifies if the plot should be framed. Default is TRUE, as in
    \code{\link{plot.default}}
  }
  \item{postscript}{
    If TRUE, (a) postscript graph(s) will be created. Otherwise, and if
    the pdf and png arguments are also FALSE, graphics will be
    displayed with a screen graphics device.
    }
  \item{pdf}{
    If TRUE, (a) pdf graph(s) will be created. Otherwise, and if
    the postscript, and png arguments are also FALSE, graphics will be
    displayed with a screen graphics device.
    }
  \item{png}{
    If TRUE, (a) png graph(s) will be created. Otherwise, and if the
    postscript and pdf arguments are also FALSE, graphics will be displayed
    with a screen graphics device.
    }
  \item{bw}{
    A boolean deciding if the plots will be black and white or not. Default
    is TRUE.
    }
  \item{pointsize}{
    The pointsize used for pdf, png and postscript graphics.
    }
  \item{colors}{
    This is a vector of colors, defaulting to 1:8, used for plotting the data.
    }
  \item{ltys}{
    This is a vector of line types for the dose-response models, defaulting to 1:8.
    }
  \item{pchs}{
    This is a vector of character types for the data. The default uses built-in
    characters 1:n with n being the number of substances for which data are plotted
    for overlays, or always 1 when no overlay plot is generated.
    }
  \item{lpos}{
    An optional argument defaulting to "topright" specifying the position
    of the legend by being passed to the legend function. See the help for the
    legend function for all possiblities.}
  \item{devoff}{
    If set to FALSE, the device will not be closed after creation of an overlay
    pdf, png or postscript graph, so texts and other elements can
    be added.
    }
}
\value{
  \item{results}{
    You will get plots of data and/or the fitted dose-response curves, on the
    screen and/or as postscript/pdf/png files, depending on the parameters.
    }
}
\examples{
data(antifoul)
r <- drfit(antifoul)
drplot(r,antifoul)
}
\note{There is a demo for each dataset that can be accessed by
  \code{demo(dataset)}}
\author{
  Johannes Ranke
  \email{jranke@uni-bremen.de}
}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
