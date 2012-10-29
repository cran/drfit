pkgname <- "drfit"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('drfit')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("IM1xIPC81")
### * IM1xIPC81

flush(stderr()); flush(stdout())

### Name: IM1xIPC81
### Title: Dose-Response data for 1-methyl-3-alkylimidazolium
###   tetrafluoroborates in IPC-81 cells
### Aliases: IM1xIPC81
### Keywords: datasets

### ** Examples

  ## Not run: demo(IM1xIPC81)



cleanEx()
nameEx("IM1xVibrio")
### * IM1xVibrio

flush(stderr()); flush(stdout())

### Name: IM1xVibrio
### Title: Dose-Response data for 1-methyl-3-alkylimidazolium
###   tetrafluoroborates in V. fischeri
### Aliases: IM1xVibrio
### Keywords: datasets

### ** Examples

  ## Not run: demo(IM1xVibrio)



cleanEx()
nameEx("XY")
### * XY

flush(stderr()); flush(stdout())

### Name: XY
### Title: Dose-Response data for two substances X and Y
### Aliases: XY
### Keywords: datasets

### ** Examples

  ## Not run: demo(XY)



cleanEx()
nameEx("antifoul")
### * antifoul

flush(stderr()); flush(stdout())

### Name: antifoul
### Title: Dose-Response data for TBT and Zink Pyrithione in IPC-81 cells
### Aliases: antifoul
### Keywords: datasets

### ** Examples

  ## Not run: demo(antifoul)



cleanEx()
nameEx("checkexperiment")
### * checkexperiment

flush(stderr()); flush(stdout())

### Name: checkexperiment
### Title: Check raw data from a specified experiment or microtiter plate
### Aliases: checkexperiment checkplate
### Keywords: database

### ** Examples

# Check plate number 3 in the cytotox database
## Not run: checkplate(3)



cleanEx()
nameEx("checksubstance")
### * checksubstance

flush(stderr()); flush(stdout())

### Name: checksubstance
### Title: Check raw data for a specified substance
### Aliases: checksubstance
### Keywords: database internal

### ** Examples

# Check substance IM14 BF4 in the cytotox database
## Not run: checksubstance("IM14 BF4")



cleanEx()
nameEx("drdata")
### * drdata

flush(stderr()); flush(stdout())

### Name: drdata
### Title: Get dose-response data via RODBC
### Aliases: drdata
### Keywords: IO database

### ** Examples

# Get cytotoxicity data for Tributyltin and zinc pyrithione, tested with IPC-81
# cells
## Not run: drdata(c("TBT","ZnPT2"))



cleanEx()
nameEx("drfit-package")
### * drfit-package

flush(stderr()); flush(stdout())

### Name: drfit-package
### Title: Dose-response data evaluation
### Aliases: drfit-package
### Keywords: package models regression nonlinear

### ** Examples

data(antifoul)
r <- drfit(antifoul)
format(r,digits=2)
drplot(r,antifoul,overlay=TRUE,bw=FALSE)



cleanEx()
nameEx("drfit")
### * drfit

flush(stderr()); flush(stdout())

### Name: drfit
### Title: Fit dose-response models
### Aliases: drfit
### Keywords: models regression nonlinear

### ** Examples

data(antifoul)
r <- drfit(antifoul)
format(r,digits=2)



cleanEx()
nameEx("drplot")
### * drplot

flush(stderr()); flush(stdout())

### Name: drplot
### Title: Plot dose-response models
### Aliases: drplot
### Keywords: models regression nonlinear

### ** Examples

data(antifoul)
r <- drfit(antifoul)
## Not run: drplot(r,antifoul)



cleanEx()
nameEx("pyrithione")
### * pyrithione

flush(stderr()); flush(stdout())

### Name: pyrithione
### Title: Cytotoxicity data for different pyrithionates and related
###   species
### Aliases: pyrithione
### Keywords: datasets

### ** Examples

  ## Not run: demo(pyrithione)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
