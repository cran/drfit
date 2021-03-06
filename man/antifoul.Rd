\name{antifoul}
\docType{data}
\alias{antifoul}
\title{Dose-Response data for TBT and Zink Pyrithione in IPC-81 cells}
\description{
  This data set shows the response of the rat leukaemic cell line IPC-81 to
  dilution series of tributyltin chloride (TBT) and Zink Pyrithione as retrieved
  from the "cytotox" database of the UFT Department of Bioorganic Chemistry on
  February 25, 2004
}
\usage{data(antifoul)}
\format{
  A dataframe containing 135 and 81 data points for concentrations and responses
  for TBT and Zink Pyrithione, respectively. Some additional columns from the
  database are also present.
}
\examples{
rantifoul.ED50 <- drfit(antifoul,
                        linlogit = TRUE, logit = TRUE, weibull = TRUE,
                        chooseone = FALSE,
                        showED50 = TRUE, EDx = c(10))
print(rantifoul.ED50, digits = 5)

rantifoul.drc <- drcfit(antifoul,
                        linlogit = TRUE, logit = TRUE, weibull = TRUE,
                        chooseone = FALSE,
                        showED50 = TRUE, EDx = c(10))
print(rantifoul.drc, digits = 5)
}

\keyword{datasets}
