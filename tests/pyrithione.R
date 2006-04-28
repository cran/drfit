library(drfit)
data(pyrithione)
rpyr <- drfit(pyrithione,linlogit=TRUE,linlogitWrong=c("MSPT","MSPHI"))
rpyr
