\name{checkexperiment}
\alias{checkexperiment}
\alias{checkplate}
\title{Check raw data from a specified experiment or microtiter plate}
\description{
  Report metadata from a specified experiment or microtiter plate from a
  specified database, box plot controls, and plot the dose-response data.
}
\usage{
  checkplate(id, db = "cytotox")
  checkexperiment(id, db = "ecotox", endpoint = "\%")
}
\arguments{
  \item{id}{
    The id of the experiment or the plate identifying it within the database.}
  \item{db}{
    The database to be used. Currently, the microtiter plate databases
    "cytotox", "enzymes" of the UFT Department of Bioorganic Chemistry are
    supported, as well as the database of ecotoxicity experiments "ecotox".}
  \item{endpoint}{
    The endpoint that should be used for selecting the data. Only important if
    the database "ecotox" is used. Defaults to "\%".}
}
\value{
  The function lists a report and shows two graphs.
}
\examples{
# Check plate number 3 in the cytotox database
\dontrun{checkplate(3)}
}
\author{
  Johannes Ranke
  \email{jranke@uni-bremen.de}
  \url{http://www.uft.uni-bremen.de/chemie/ranke}
}
\keyword{database}
