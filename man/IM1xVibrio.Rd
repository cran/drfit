\name{IM1xVibrio}
\docType{data}
\alias{IM1xVibrio}
\title{Dose-Response data for 1-methyl-3-alkylimidazolium tetrafluoroborates in V. fischeri}
\description{
  This is the raw data documenting the influence of the alkyl chain length in 3
  position on the toxicity to the marine luminescent bacteria \emph{V.
  fischeri}. The substances are named according to the UFT naming scheme of
  these substances. 
  IM13 BF4 means 1-methyl-3-propylimidazolium tetrafluoroborate, 
  IM14 BF4 means 1-methyl-3-butylimidazolium tetrafluoroborate and 
  IM1-10 BF4 means 1-methyl-3-decylimidazolium tetrafluoroborate.
}
\usage{data(IM1xVibrio)}
\format{
  A dataframe containing the data as required for the \code{\link{drfit}}
  function. Additional columns contain the species tested (luminescent bacteria
  Vibrio fischeri, \code{organism}), and a field specifying if the data is
  regarded valid (\code{ok}).
}
\examples{
  \dontrun{demo(IM1xVibrio)}
}
\source{
  Ranke J, Mölter K, Stock F, Bottin-Weber U, Poczobutt J, Hoffmann J,
  Ondruschka B, Filser J, Jastorff B (2004) Biological effects of imidazolium
  ionic liquids with varying chain lenghts in acute Vibrio fischeri and WST-1
  cell viability assays. Ecotoxicology and Environmental Safety 58(3) 396-404
}
\keyword{datasets}
