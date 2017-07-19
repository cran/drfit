\name{drdata}
\alias{drdata}
\title{Get dose-response data via RODBC}
\description{
  Get dose-response data from an adequate ODBC data source
}
\usage{
  drdata(substances, experimentator = "\%", db = "cytotox", celltype = "IPC-81",
    enzymetype="AChE", organism="Vibrio fischeri", endpoint="Luminescence",
    whereClause = "1", ok = "'ok','no fit'")
}
\arguments{
  \item{substances}{
    A string or an array of strings with the substance names for
    which dose-response data is to be retrieved. Matching is case-sensitive.}
  \item{experimentator}{
    The name of the experimentator whose data is to be used. Default is "%",
    which means that data from all experimentators are retrieved.}
  \item{db}{
    The database to be used. Currently, the databases "cytotox", "enzymes"
    and "ecotox" of the UFT Department of Bioorganic Chemistry are
    supported (default is "cytotox").}
  \item{celltype}{
    Currently, only data for IPC-81, C6, NB4, HeLa, Jurkat and U937 are
    supported.}
  \item{enzymetype}{
    Currently, only data for AChE, GR and GST are supported.}
  \item{organism}{
    The organism that was exposed to the chemical. Only important if the database
    "ecotox" is used. Defaults to "Vibrio fischeri".}
  \item{endpoint}{
    The endpoint that should be used for selecting the data. Only important if
    the database "ecotox" is used. Defaults to "Luminescence".}
  \item{whereClause}{
    With this argument, additional conditions for the SQL query can be set,
    e.g. "plate != 710" (i.e. "Do not retrieve data for plate 710"). The
    default is 1 (in SQL syntax this means TRUE).}
  \item{ok}{
    With the default value "'ok','no fit'", only data that has been checked and
    set to "ok" or "no fit" in the database is retrieved. The argument "no fit"
    will result in not using the data for fitting, but it will be plotted.
    Another sensible argument would be "'ok','no fit','?'", in order to
    additionally retrieve data which has not yet been checked.}
}
\value{
  \item{data}{
    A data frame with a factor describing the dose levels, the numeric dose
    levels and a numeric column describing the response, as well as the entries
    for plate, experimentator, performed (date of test performance), celltype,
    unit (of the dose/concentration), and for the ok field in the database.}
}
\details{
  The function is currently only used for retrieving data from the
  mysql database "cytotox" of the UFT Department of Bioorganic Chemistry.
  Access to this database is limited to UFT staff.  Additionally to the
  installation of the RODBC package, it is required to set up a ODBC data
  source with the name "cytotox", using an ODBC driver for mysql, probably
  myODBC. Then, under Unix, you can use iodbc or unixodbc for setting up the
  respective data source with data source name (DSN) "cytotox". For my
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
  }
}
\examples{
# Get cytotoxicity data for Tributyltin and zinc pyrithione, tested with IPC-81
# cells
\dontrun{drdata(c("TBT","ZnPT2"))}
}
\author{
  Johannes Ranke
  \email{jranke@uni-bremen.de}
  \url{http://www.uft.uni-bremen.de/chemie/ranke}
}
\keyword{IO}
\keyword{database}
