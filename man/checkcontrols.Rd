\name{checkcontrols}
\alias{checkcontrols}
\title{Check raw data from range of plates or experiments}
\description{
  Report control data from a specified range of microtiter plates or
  experiments from a specified database.
}
\usage{
  checkcontrols(last = 10, id = NULL, db = "cytotox", 
                celltype = "IPC-81", enzymetype = "AChE", 
                organism = "Vibrio fischeri",
                endpoint = "\%", qcc = c("R", "xbar"))
}
\arguments{
  \item{last}{
    If \code{id} is not specified, this argument specifies the number of
    plates/experiments that should be evaluated. The plates/experiments with
    the hightest numbers are selected.
  }
  \item{id}{
    A numeric vector of integers, specifying the plates or experiments explicitly.
  }
  \item{db}{
    The database to be used. Currently, the microtiter plate databases
    "cytotox", "enzymes" of the UFT Department of Bioorganic Chemistry are
    supported, as well as the database of ecotoxicity experiments "ecotox".
  }
  \item{celltype}{
    Only important if database "cytotox" is used.  Data for
    "IPC-81", "C6", "NB4", "HeLa", "Jurkat" and "U937" are available.
  }
  \item{enzymetype}{
    Only important if database "enzymes" is used.
    Data for "AChE", "GR" and "GST" are available.
  }
  \item{organism}{
    Only important if database "ecotox" is used. 
  }
  \item{endpoint}{
    The endpoint that should be used for selecting the data. Only important if
    the database "ecotox" is used. Defaults to "\%".}
  \item{qcc}{
    The type of quality control charts to be plotted. By default, an R chart
    (showing ranges of control values within plates/experiments), and an
    xbar chart (showing means) are generated.
  }
}
\value{
  The function lists a report and shows two graphs.
}
\examples{
\dontrun{checkcontrols(15)}
}
\author{
  Johannes Ranke
  \email{jranke@uni-bremen.de}
  \url{http://www.uft.uni-bremen.de/chemie/ranke}
}
\keyword{database}
