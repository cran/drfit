drdata <- function(substances, experimentator = "%", db = "cytotox",
    celltype="IPC-81",enzymetype="AChE",
    organism="Vibrio fischeri",endpoint="Luminescence",whereClause="1",
    ok="'ok','no fit'")
{
    library(RODBC) 
    channel <- odbcConnect(db,uid="cytotox",pwd="cytotox",case="tolower")
    slist <- paste(substances,collapse="','")
    if (db == "cytotox") {
        responsetype <- "viability"
        testtype <- "celltype"
        type <- celltype
    } else {
        if (db == "enzymes") {
            responsetype <- "activity"
            testtype <- "enzyme"
            type <- enzymetype
        } else {
            responsetype <- "response"
            testtype <- "organism"
            type <- organism
        }
    }
        
    query <- paste("SELECT conc,",responsetype,",unit,experimentator,substance,",testtype,
        ",ok FROM ", db, " WHERE substance IN ('",
        slist,"') AND experimentator LIKE '",
        experimentator,"' AND ",testtype," LIKE '",
        type,"' AND ",
        whereClause," AND ok in (",
        ok,")",sep="")
    if (db == "ecotox") query <- paste(query," AND type LIKE '",endpoint,"'",sep="")
    data <- sqlQuery(channel,query)
    odbcClose(channel)
    names(data)[[1]] <- "dose"
    names(data)[[2]] <- "response"
    data$dosefactor <- factor(data$dose)
    data$substance <- factor(data$substance,levels=substances)
    return(data)
}
    
linlogitf <- function(x,k,f,mu,b)
{
    k*(1 + f*x) / (1 + ((2*f*(10^mu) + 1) * ((x/(10^mu))^b)))
}

drfit <- function(data, startlogED50 = NA, chooseone=TRUE,
        probit = TRUE, logit = FALSE, weibull = FALSE,
        linlogit = FALSE, linlogitWrong = NA, allWrong = NA,
        s0 = 0.5, b0 = 2, f0 = 0)
{
    if(!is.null(data$ok)) data <- subset(data,ok!="no fit") # Don't use data where ok 
                                                            # was set to "no fit"
    substances <- levels(data$substance)

    ri <- rix <- 0                        # ri is the index over the result rows
                                          # rix is used later to check if any
                                          # model result was appended
    rsubstance <- array()                 # the substance names in the results
    rndl <- vector()                      # number of dose levels
    rn <- vector()                        # mean number of replicates in each dose level
    runit <- vector()                     # vector of units for each result row
    rlhd <- rlld <- vector()              # highest and lowest doses tested
    mtype <- array()                      # the modeltypes
    sigma <- array()                      # the standard deviation of the residuals
    logED50 <- vector()
    stderrlogED50 <- vector()
    a <- b <- c <- vector()

    splitted <- split(data,data$substance)
    for (i in substances) {
        tmp <- splitted[[i]]
        fit <- FALSE
        if (length(tmp) != 0) {
            unit <- levels(as.factor(as.vector(tmp$unit)))
            cat("\n",i,": Fitting data...\n",sep="")
        } else {
            unit <- ""
            cat("\n",i,": No data\n",sep="")
        }
        if (length(unit) > 1) {
            cat("More than one unit for substance ",i,", halting\n\n",sep="")
            break
        }
        if (length(tmp$response) == 0) {
            nodata = TRUE
        } else {
            nodata = FALSE
        }
        rix <- ri
        if (nodata) {
            n <- ndl <- 0
        } else {
            ndl <- length(levels(factor(tmp$dose)))
            n <- round(length(tmp$response)/ndl)
            if (is.na(startlogED50[i])){
                w <- 1/abs(tmp$response - 0.3)
                startlogED50[[i]] <- sum(w * log10(tmp$dose))/sum(w)
            }
            highestdose <- max(tmp$dose)
            lowestdose <- min(tmp$dose)
            lhd <- log10(highestdose)
            lld <- log10(lowestdose)
            responseathighestdose <- mean(subset(tmp,dose==highestdose)$response)
            if (responseathighestdose < 0.5) {
                inactive <- FALSE

                if (linlogit && 
                    length(subset(linlogitWrong,linlogitWrong == i))==0 &&
                    length(subset(allWrong,allWrong == i))==0) {
                    m <- try(nls(response ~ linlogitf(dose,1,f,logED50,b),
                            data=tmp,
                            start=list(f=f0,logED50=startlogED50[[i]],b=b0)))
                    if (!inherits(m, "try-error")) {
                        fit <- TRUE
                        ri <- ri + 1
                        s <- summary(m)
                        sigma[[ri]] <- s$sigma
                        rsubstance[[ri]] <- i
                        rndl[[ri]] <- ndl
                        rn[[ri]] <- n
                        runit[[ri]] <- unit
                        rlld[[ri]] <- log10(lowestdose)
                        rlhd[[ri]] <- log10(highestdose)
                        logED50[[ri]] <- coef(m)[["logED50"]]
                        if (logED50[[ri]] > rlhd[[ri]]) {
                            mtype[[ri]] <- "no fit"
                            logED50[[ri]] <- NA
                            stderrlogED50[[ri]] <- NA
                            a[[ri]] <- NA
                            b[[ri]] <- NA
                            c[[ri]] <- NA
                        } else {
                            mtype[[ri]] <- "linlogit"
                            stderrlogED50[[ri]] <- s$parameters["logED50","Std. Error"]
                            a[[ri]] <- coef(m)[["logED50"]]
                            b[[ri]] <- coef(m)[["b"]]
                            c[[ri]] <- coef(m)[["f"]]
                        }
                    }
                }

                if (probit &&
                    length(subset(allWrong,allWrong == i))==0) {
                    m <- try(nls(response ~ pnorm(-log10(dose),-logED50,scale),
                                data=tmp,
                                start=list(logED50=startlogED50[[i]],scale=1)))
                    if (chooseone==FALSE || fit==FALSE) {
                        if (!inherits(m, "try-error")) {
                            fit <- TRUE
                            ri <- ri + 1
                            s <- summary(m)
                            sigma[[ri]] <- s$sigma
                            rsubstance[[ri]] <- i
                            rndl[[ri]] <- ndl
                            rn[[ri]] <- n
                            runit[[ri]] <- unit
                            rlld[[ri]] <- log10(lowestdose)
                            rlhd[[ri]] <- log10(highestdose)
                            logED50[[ri]] <- coef(m)[["logED50"]]
                            c[[ri]] <- NA
                            if (logED50[[ri]] > rlhd[[ri]]) {
                                mtype[[ri]] <- "no fit"
                                logED50[[ri]] <- NA
                                stderrlogED50[[ri]] <- NA
                                a[[ri]] <- NA
                                b[[ri]] <- NA
                            } else {
                                mtype[[ri]] <- "probit"
                                stderrlogED50[[ri]] <- s$parameters["logED50","Std. Error"]
                                a[[ri]] <- coef(m)[["logED50"]]
                                b[[ri]] <- coef(m)[["scale"]]
                            }
                        }
                    }
                }

                if (logit &&
                    length(subset(allWrong,allWrong == i))==0) {
                    m <- try(nls(response ~ plogis(-log10(dose),-logED50,scale),
                            data=tmp,
                            start=list(logED50=startlogED50[[i]],scale=1)))
                    if (chooseone==FALSE || fit==FALSE) {
                        if (!inherits(m, "try-error")) {
                            fit <- TRUE
                            ri <- ri + 1
                            s <- summary(m)
                            sigma[[ri]] <- s$sigma
                            rsubstance[[ri]] <- i
                            rndl[[ri]] <- ndl
                            rn[[ri]] <- n
                            runit[[ri]] <- unit
                            rlld[[ri]] <- log10(lowestdose)
                            rlhd[[ri]] <- log10(highestdose)
                            logED50[[ri]] <- a[[ri]] <- coef(m)[["logED50"]]
                            b[[ri]] <- coef(m)[["scale"]]
                            c[[ri]] <- NA
                            if (logED50[[ri]] > rlhd[[ri]]) {
                                mtype[[ri]] <- "no fit"
                                logED50[[ri]] <- NA
                                stderrlogED50[[ri]] <- NA
                                a[[ri]] <- NA
                                b[[ri]] <- NA
                            } else {
                                mtype[[ri]] <- "logit"
                                stderrlogED50[[ri]] <- s$parameters["logED50","Std. Error"]
                            }
                        }
                    }
                }

                if (weibull &&
                    length(subset(allWrong,allWrong == i))==0) {
                    m <- try(nls(response ~ pweibull(-log10(dose)+location,shape),
                            data=tmp,
                            start=list(location=startlogED50[[i]],shape=s0)))
                    if (chooseone==FALSE || fit==FALSE) {
                        if (!inherits(m, "try-error")) {
                            fit <- TRUE
                            ri <- ri + 1
                            s <- summary(m)
                            sigma[[ri]] <- s$sigma
                            rsubstance[[ri]] <- i
                            rndl[[ri]] <- ndl
                            rn[[ri]] <- n
                            runit[[ri]] <- unit
                            rlld[[ri]] <- log10(lowestdose)
                            rlhd[[ri]] <- log10(highestdose)
                            a[[ri]] <- coef(m)[["location"]]
                            b[[ri]] <- coef(m)[["shape"]]
                            sqrdev <- function(logdose) {
                                (0.5 - pweibull( - logdose + a[[ri]], b[[ri]]))^2 
                            }
                            logED50[[ri]] <- nlm(sqrdev,startlogED50[[i]])$estimate
                            c[[ri]] <- NA
                            if (logED50[[ri]] > rlhd[[ri]]) {
                                mtype[[ri]] <- "no fit"
                                logED50[[ri]] <- NA
                                stderrlogED50[[ri]] <- NA
                                a[[ri]] <- NA
                                b[[ri]] <- NA
                            } else {
                                mtype[[ri]] <- "weibull"
                                stderrlogED50[[ri]] <- NA
                            }
                        }
                    }
                }


            } else {
                inactive <- TRUE
            }
        }
        if (ri == rix) {          # if no entry was appended for this substance
            ri <- ri + 1
            rsubstance[[ri]] <- i
            rndl[[ri]] <- ndl
            rn[[ri]] <- n
            if (nodata) {
                rlld[[ri]] <- rlhd[[i]] <- NA
                mtype[[ri]] <- "no data"
                runit[[ri]] <- NA
            } else {
                rlld[[ri]] <- log10(lowestdose)
                rlhd[[i]] <- log10(highestdose)
                runit[[ri]] <- unit
                if (inactive) {
                    mtype[[ri]] <- "inactive"
                } else {
                    mtype[[ri]] <- "no fit"
                }
            }
            sigma[[ri]] <- NA
            logED50[[ri]] <- NA
            stderrlogED50[[ri]] <- NA
            a[[ri]] <- NA
            b[[ri]] <- NA
            c[[ri]] <- NA
        }
    }
    results <- data.frame(rsubstance, rndl, rn, rlld, rlhd, mtype, logED50, stderrlogED50, runit, sigma, a, b)
    names(results) <- c("Substance","ndl","n","lld","lhd","mtype","logED50","std","unit","sigma","a","b")
    if (linlogit) {
        results$c <- c
    }
    return(results)
}

drplot <- function(drresults, data, dtype = "std", alpha = 0.95,
        path = "./", fileprefix = "drplot", overlay = FALSE,
        postscript = FALSE, png = FALSE, bw = TRUE,
        pointsize = 12,
        colors = 1:8, devoff=T, lpos="topright")
{
    # Check if all data have the same unit
    unitlevels <- levels(as.factor(drresults$unit))
    if (length(unitlevels) == 1) {
        unit <- unitlevels
    } else {
        unit <- "different units"
    }

    # Get the plot limits on the x-axis (log of the dose)
    if(is.data.frame(data)) {
        if (min(data$dose) == 0) {
            cat("At least one of the dose levels is 0 - this is not a valid dose.")
        } else {
            lld <- log10(min(data$dose))
        }
        lhd <- log10(max(data$dose))
        hr <- max(data$response)
        dsubstances <- levels(data$substance)    
    } else {
        lld <- min(drresults[["logED50"]],na.rm=TRUE) - 2
        lhd <- max(drresults[["logED50"]],na.rm=TRUE) + 2
        if (length(subset(drresults,mtype=="linlogit")$Substance) != 0) {
            hr <- 1.8 
        } else {
            hr <- 1.0
        }
    }

    # Prepare overlay plot if requested
    if (overlay)
    {
        if (postscript) {
            filename = paste(path,fileprefix,".eps",sep="")
            postscript(file=filename,
                    paper="special",width=7,height=7,horizontal=FALSE, pointsize=pointsize)
            cat("Created File: ",filename,"\n")
        } 
        if (png) {
            filename = paste(path,fileprefix,".png",sep="")
            png(filename=filename,
                width=500, height=500, pointsize=pointsize)
            cat("Created File: ",filename,"\n")
        }
        if (!postscript && !png) {
            get(getOption("device"))(width=7,height=7)
        }
            
        plot(0,type="n",
            xlim=c(lld - 0.5, lhd + 1),
            ylim= c(-0.1, hr + 0.2),
            xlab=paste("Decadic Logarithm of the dose in ", unit),    
            ylab="Normalized response")
    }

    # Plot the data either as raw data or as error bars
    if(is.data.frame(data)) {
        splitted <- split(data,data$substance)
        # n is the index for the dose-response curves
        n <- 0
        if (bw) colors <- rep("black",length(dsubstances))
        # Loop over the substances in the data
        for (i in dsubstances) {
            n <- n + 1
            tmp <- splitted[[i]]
            if (length(tmp$response) != 0) {
                color <- colors[[n]]
                # Prepare the single graphs if an overlay is not requested
                if (!overlay)
                {
                    if (postscript) {
                        filename = paste(path,fileprefix,sub(" ","_",i),".eps",sep="")
                        postscript(file=filename,
                                paper="special",width=7,height=7,horizontal=FALSE,pointsize=pointsize)
                        cat("Created File: ",filename,"\n")
                    } 
                    if (png) {
                        filename = paste(path,fileprefix,sub(" ","_",i),".png",sep="")
                        png(filename=filename,
                            width=500, height=500, pointsize=pointsize)
                        cat("Created File: ",filename,"\n")
                    }
                    if (!postscript && !png) {
                         get(getOption("device"))(width=7,height=7)
                    }
                        
                    plot(0,type="n",
                        xlim=c(lld - 0.5, lhd + 2),
                        ylim= c(-0.1, hr + 0.2),
                        xlab=paste("Decadic Logarithm of the dose in ", unit),    
                        ylab="Normalized response")
                }
                if (!overlay) legend(lpos, i,lty = 1, col = color, inset=0.05)
                tmp$dosefactor <- factor(tmp$dose)  # necessary because the old
                                                    # factor has all levels, not 
                                                    # only the ones tested with
                                                    # this substance
                # Plot the data, if requested
                if (dtype != "none") {
                    if (dtype == "raw") {
                        points(log10(tmp$dose),tmp$response,col=color)
                    } else {
                        splitresponses <- split(tmp$response,tmp$dosefactor)
                        means <- sapply(splitresponses,mean)
                        lengths <- sapply(splitresponses,length)
                        vars <- sapply(splitresponses,var)
                        standarddeviations <- sqrt(vars)
                    }
                    if (dtype == "std")
                    {
                        tops <- means + standarddeviations
                        bottoms <- means - standarddeviations
                    }
                    if (dtype == "conf")
                    {
                        confidencedeltas <- qt((1 + alpha)/2, lengths - 1) * sqrt(vars) 
                        tops <- means + confidencedeltas
                        bottoms <- means - confidencedeltas
                    }
                    if (dtype != "raw")
                    {
                        x <- log10(as.numeric(levels(tmp$dosefactor)))
                        segments(x,bottoms,x,tops,col=color)
                        points(x,means,col=color)
                        smidge <- 0.05
                        segments(x - smidge,bottoms,x + smidge,bottoms,col=color)
                        segments(x - smidge,tops,x + smidge,tops,col=color)
                    }
                }

                # Plot the fits, if there are any
                fits <- subset(drresults,Substance == i)
                nf <-  length(fits$Substance)  # number of fits to plot
                if (nf > 0) {
                    for (j in 1:nf)
                    {
                        logED50 <- fits[j,"logED50"]
                        mtype <- as.character(fits[j, "mtype"])
                        if (mtype == "probit") {
                            scale <- fits[j,"b"]
                            plot(function(x) pnorm(-x,-logED50,scale),lld - 0.5, lhd + 2, add=TRUE,col=color)
                        }
                        if (mtype == "logit") {
                            scale <- fits[j,"b"]
                            plot(function(x) plogis(-x,-logED50,scale),lld - 0.5, lhd + 2, add=TRUE,col=color)
                        }
                        if (mtype == "weibull") {
                            location <- fits[j,"a"]
                            shape <- fits[j,"b"]
                            plot(function(x) pweibull(-x+location,shape),lld - 0.5, lhd + 2, add=TRUE,col=color)
                        }
                        if (mtype == "linlogit") {
                            plot(function(x) linlogitf(10^x,1,fits[j,"c"],fits[j,"logED50"],fits[j,"b"]),
                                lld - 0.5, lhd + 2,
                                add=TRUE,col=color)
                        }
                    }
                }
                if (!overlay && (postscript || png)) dev.off()
            } else {
                cat("No data for ",i,"\n")
            }
        }
    }
    if (overlay) legend(lpos, dsubstances,lty = 1, col = colors, inset=0.05)
    if (overlay && (postscript || png)) {
        if (devoff) {
            dev.off()
        }
    }
}

checkplate <- function(plate,db="cytotox")
{
    library(RODBC) 
    channel <- odbcConnect(db,uid="cytotox",pwd="cytotox",case="tolower")

    if (db == "cytotox") {
        responsetype <- "viability"
        testtype <- "celltype"
    } else {
        responsetype <- "activity"
        testtype <- "enzyme"
    }
        
    platequery <- paste("SELECT experimentator,substance,",testtype,",conc,unit,",responsetype,",performed,ok",
        "FROM ",db," WHERE plate=", plate)

    controlquery <- paste("SELECT type,response FROM controls WHERE plate=",plate)
    
    platedata <- sqlQuery(channel,platequery)
    controldata <- sqlQuery(channel,controlquery)

    odbcClose(channel)

    if (length(platedata$experimentator) < 1) {
        cat("There is no response data for plate ",plate," in database ",db,"\n")
    } else {
        platedata$experimentator <- factor(platedata$experimentator)
        platedata$type <- factor(platedata[[testtype]])
        platedata$substance <- factor(platedata$substance)
        platedata$unit <- factor(platedata$unit)
        platedata$performed <- factor(platedata$performed)
        platedata$ok <- factor(platedata$ok)
        
        blinds <- subset(controldata,type=="blind")
        controls <- subset(controldata,type=="control")
        
        numberOfBlinds <- length(blinds$response)
        numberOfControls <- length(controls$response)
        meanOfBlinds <- mean(blinds$response)
        meanOfControls <- mean(controls$response)
        stdOfBlinds <- sd(blinds$response)
        stdOfControls <- sd(controls$response)
        percentstdOfcontrols <-stdOfControls *100/meanOfControls
        
        cat("Plate ",plate," from database ",db,"\n",
            "\tExperimentator: ",levels(platedata$experimentator),"\n",
            "\tType(s): ",levels(platedata$type),"\n",
            "\tPerformed on : ",levels(platedata$performed),"\n",
            "\tSubstance(s): ",levels(platedata$substance),"\n",
            "\tConcentration unit: ",levels(platedata$unit),"\n",
            "\tOK: ",levels(platedata$ok),"\n",
            "\t\tNumber \tMean \t\tStandard Deviation \t% Standard Deviation \n",
            "\tblind\t",numberOfBlinds,"\t",meanOfBlinds,"\t",stdOfBlinds,"\n",
            "\tcontrol\t",numberOfControls,"\t",meanOfControls,"\t",stdOfControls,"\t\t",percentstdOfcontrols,"\n")
        
        par(ask=TRUE)
        
        boxplot(blinds$response,controls$response,names=c("blinds","controls"),ylab="Response",main=paste("Plate ",plate))
        
        drdata <- platedata[c(2,4,6)]
        drdata$substance <- factor(drdata$substance)
        substances <- levels(drdata$substance)
       
        plot(log10(drdata$conc),drdata$viability,
            xlim=c(-2.5, 4.5), 
            ylim= c(-0.1, 2), 
            xlab=paste("decadic logarithm of the concentration in ",levels(platedata$unit)),
            ylab=responsetype)
        
        drdatalist <- split(drdata,drdata$substance)
        
        for (i in 1:length(drdatalist)) {
            points(log10(drdatalist[[i]]$conc),drdatalist[[i]][[responsetype]],col=i);
        }

        legend("topleft",substances, pch=1, col=1:length(substances), inset=0.05)
        title(main=paste("Plate ",plate," - ",levels(platedata$experimentator)," - ",levels(platedata$type)))
    }
}

checksubstance <- function(substance,db="cytotox",experimentator="%",celltype="%",enzymetype="%",whereClause="1",ok="%") 
{
    library(RODBC) 
    channel <- odbcConnect(db,uid="cytotox",pwd="cytotox",case="tolower")

    if (db == "cytotox") {
        responsetype <- "viability"
        testtype <- "celltype"
        type <- celltype
    } else {
        responsetype <- "activity"
        testtype <- "enzyme"
        type <- enzymetype
    }
    query <- paste("SELECT experimentator,substance,",testtype,",plate,conc,unit,",responsetype,",ok",
        " FROM ",db," WHERE substance LIKE '",
        substance,"' AND experimentator LIKE '",
        experimentator,"' AND ",testtype," LIKE '",
        type,"' AND ",
        whereClause," AND ok LIKE '",ok,"'",sep="")

    data <- sqlQuery(channel,query)
    odbcClose(channel)
    
    data$experimentator <- factor(data$experimentator)    
    data$substance <- factor(data$substance)
    substances <- levels(data$substance)    
    data$type <- factor(data[[testtype]])                
    data$plate <- factor(data$plate)                        
    plates <- levels(data$plate)
    concentrations <- split(data$conc,data$conc)
    concentrations <- as.numeric(names(concentrations))
    data$unit <- factor(data$unit)                        
    data$ok <- factor(data$ok)
    
    if (length(plates)>6) {
        palette(rainbow(length(plates)))      
    }
 
    plot(log10(data$conc),data[[responsetype]],
        xlim=c(-2.5, 4.5),                                                                  
        ylim= c(-0.1, 2),                                                                 
        xlab=paste("decadic logarithm of the concentration in ",levels(data$unit)),    
        ylab=responsetype)  
        
    platelist <- split(data,data$plate)
   
    for (i in 1:length(platelist)) {    
        points(log10(platelist[[i]]$conc),platelist[[i]][[responsetype]],col=i);          
    }       
    
    legend("topleft", plates, pch=1, col=1:length(plates), inset=0.05)
    title(main=paste(substance," - ",levels(data$experimentator)," - ",levels(data$type)))
 
    cat("Substanz ",substance,"\n",
        "\tExperimentator(s):",levels(data$experimentator),"\n",
        "\tType(s):\t",levels(data$type),"\n",
        "\tSubstance(s):\t",levels(data$substance),"\n",
        "\tPlate(s):\t",plates,"\n\n")
}
