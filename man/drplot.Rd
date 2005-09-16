\name{drplot}
\alias{drplot}
\title{Plot dose-response models}
\description{
	Produce graphics of dose-response data and dose-response relationships 
  either combined or separately, for one or more substances.
}
\usage{
  drplot(drresults, data, dtype, alpha, path, fileprefix, overlay,
    postscript, png, bw, pointsize, colors, devoff, lpos)
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
  \item{path}{
    The path where graphic files should be put if any are produced. Defaults
    to "./" i.e. the current working directory of R.
    }
  \item{fileprefix}{
    A string which will form the beginning of each filename, if graphic files are 
    created. Defaults to "drplot".
    }
  \item{overlay}{
    If TRUE, all output will be put into one graph, otherwise a separate graph
    will be created for each substance. In the latter case, on-screen display
    (postscript=FALSE) only works correctly for data plots. Dose-response models
    will all be put into the last graph in this case.
    }
  \item{postscript}{
    If TRUE, (a) postscript graph(s) will be created. Otherwise, and if 
    the png argument is also FALSE, graphics will be
    displayed with a screen graphics device.
    }
  \item{png}{
    If TRUE, (a) png graph(s) will be created. Otherwise, and if the
    postscript argument is also FALSE, graphics will be displayed with a
    screen graphics device.
    }
  \item{bw}{
    A boolean deciding if the plots will be black and white or not. Default 
    is TRUE.
    }
  \item{pointsize}{
    The pointsize used for png and postscript graphics.
    }
  \item{colors}{
    This is a vector of colors, defaulting to 1:8, used for plotting the data.
    }
  \item{lpos}{
    An optional argument defaulting to "topright" specifying the position
    of the legend by being passed to the legend function. See the help for the
    legend function for all possiblities.} 
  \item{devoff}{
    If set to FALSE, the closing of the device after creation of an overlay png
    or postscript graph will be left out, so texts and other elements can be 
    added to the graph.
    }
}
\value{
  \item{results}{
    You will get plots of data and/or the fitted dose-response curves, on the
    screen and/or as postscript files, depending on the parameters.
    }
    
} 
\note{
  Be sure all devices are closed (e.g. by calling \code{dev.off()}) before
  calling \code{drplot} again after a failure.
}
\examples{
data(antifoul)
r <- drfit(antifoul)
\dontrun{drplot(r,antifoul)}
}
\author{
  Johannes Ranke 
  \email{jranke@uni-bremen.de} 
  \url{http://www.uft.uni-bremen.de/chemie/ranke}
}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
