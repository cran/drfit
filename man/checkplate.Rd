\name{checkplate}
\alias{checkplate}
\title{Check raw data from a specified microtiter plate}
\description{
  Report metadata from a specified microtiter plate from a specified database, box 
  plot positive and negative (blind) controls, and show the response data on the 
  plate.
}
\usage{
  checkplate(plate,db="cytotox")
}
\arguments{
  \item{plate}{
    The number of the plate identifying it within the database.}
  \item{db}{
    The database to be used. Currently, the databases "cytotox" and "enzymes"
    of the UFT Department of Bioorganic Chemistry are supported (default is
    "cytotox").} 
}
\value{
  The function lists a report and shows two graphs.
}
\examples{
# Check plate number 1 in the cytotox database
\dontrun{data <- checkplate(1)}
}
\author{
  Johannes Ranke 
  \email{jranke@uni-bremen.de} 
  \url{http://www.uft.uni-bremen.de/chemie/ranke}
}
\keyword{database}
