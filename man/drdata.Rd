\name{drdata}
\alias{drdata}
\title{Get dose-response data}
\description{
	Get dose-response data from a remote mysql server 
}
\usage{
  drdata(substances, experimentator = "\%", db = "cytotox", celltype = "IPC-81", 
    whereClause = "1", ok = "'ok'")
}
\arguments{
  \item{substances}{
    A string or an array of strings with the substance names for
    which dose-response data is to be retrieved.}
  \item{experimentator}{
    The name of the experimentator whose data is to be used.}
  \item{db}{
    The database to be used. Currently only "cytotox" of the UFT Department of
    Bioorganic Chemistry is supported.}
  \item{celltype}{
    Currently, only data for IPC-81, C6, NB4, HeLa, Jurkat and U937 are supported.}
  \item{whereClause}{
    With this argument, additional conditions for the SQL query can be set, 
    e.g. "where plate != 710". The default is 1 (in SQL syntax this means TRUE).}
  \item{ok}{
    With the default value "'ok'", only data that has been checked and set to "ok"
    in the database is retrieved. Another sensible argument would be "'ok','?'", in 
    order to additionally retrieve data which has not yet been checked.}
}
\value{
  \item{data}{
    A data frame with a factor describing the dose levels, the numeric dose levels
    and a numeric column describing the response, as well as the entries for
    plate, experimentator, performed (date of test performance), celltype, unit
    (of the dose/concentration), and for the ok field in the database.}
}
\details{
  The function is currently only used for retrieving data from the
  mysql database "cytotox" of the UFT Department of Bioorganic Chemistry.
  Additionally to the installation of the RODBC package, it is required to set
  up a ODBC data source with the name "cytotox", using an ODBC driver for mysql, 
  probably myODBC. Then, under Unix, you can use iodbc or unixodbc for setting
  up the respective data source with data source name (DSN) "cytotox". For my
  setting using unixodbc, I am using the file \file{/etc/odbcinst.ini}
  containing: 
  \tabular{lll}{
    [MySQL] \tab \tab \cr
    Description \tab = \tab MySQL driver for ODBC \cr
    Driver \tab = \tab /usr/local/lib/libmyodbc.so \cr
    Setup \tab = \tab /usr/lib/odbc/libodbcmyS.so \cr
  }
  and the file \file{/etc/odbc.ini} containing:
  \tabular{lll}{
    [cytotox] \tab \tab \cr
    Description \tab = \tab Cytotoxicity database of the department of bioorganic chemistry, UFT Bremen \cr
    Driver \tab = \tab MySQL \cr
    Trace \tab = \tab Yes  \cr
    TraceFile \tab = \tab /tmp/odbc.log  \cr
    Database \tab = \tab cytotox  \cr
    Server \tab = \tab eckehaat  \cr
    Port \tab = \tab 3306  \cr
  }.
}
\examples{
# Get cytotoxicity data for Tributyltin and zinc pyrithione, tested with IPC-81 cells
\dontrun{data <- drdata(c("TBT","Zn Pyrithion"))}
}
\author{
  Johannes Ranke 
  \email{jranke@uni-bremen.de} 
  \url{http://www.uft.uni-bremen.de/chemie/ranke}
}
\keyword{IO}
\keyword{database}
