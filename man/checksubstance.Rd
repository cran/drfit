\name{checksubstance}
\alias{checksubstance}
\title{Check raw data for a specified substance}
\description{
  Report metadata for a specified substance from a specified database, and plot
  the data.
}
\usage{
  checksubstance(substance,db="cytotox",experimentator="\%",celltype="\%",enzymetype="\%",whereClause="1",ok="\%")
}
\arguments{
  \item{substance}{
    The name of the substance identifying it within the database.}
  \item{db}{
    The database to be used. Currently, the databases "cytotox" and "enzymes"
    of the UFT Department of Bioorganic Chemistry are supported (default is
    "cytotox").} 
  \item{experimentator}{
    The name of the experimentator whose data is to be used. Default is "\%",
    which means that data from all experimentators are shown.}
  \item{celltype}{
    Currently, only data for IPC-81, C6, NB4, HeLa, Jurkat and U937 are
    supported.  Default is "\%", i.e. data for any cell type will be
    displayed.}
  \item{enzymetype}{
    Currently, only data for AChE, GR and GST are supported. The default value
    is "\%", i.e. data for any enzyme type will be displayed.}
  \item{whereClause}{
    With this argument, additional conditions for the SQL query can be set, 
    e.g. "where plate != 710". The default is 1 (in SQL syntax this means TRUE).}
  \item{ok}{
    With the default value "\%", all data in the database is retrieved for the 
    specified substance.}
}
\value{
  The function lists a report and shows one graph.
}
\examples{
# Check substance IM14 BF4 in the cytotox database
\dontrun{checksubstance("IM14 BF4")}
}
\author{
  Johannes Ranke 
  \email{jranke@uni-bremen.de} 
  \url{http://www.uft.uni-bremen.de/chemie/ranke}
}
\keyword{database}
\keyword{internal}
