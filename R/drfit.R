drdata <- function(substances, experimentator = "%", db = "cytotox",
    celltype="IPC-81",enzymetype="AChE",whereClause="1",
    ok="'ok'")
{
    library(RODBC) 
    channel <- odbcConnect(db,uid="cytotox",pwd="cytotox",case="tolower")
    slist <- paste(substances,collapse="','")
    if (db == "cytotox") {
        responsetype <- "viability"
        testtype <- "celltype"
        type <- celltype
    } else {
        responsetype <- "activity"
        testtype <- "enzyme"
        type <- enzymetype
    }
        
    query <- paste("SELECT conc,",responsetype,",unit,experimentator,substance,",testtype,
        ",plate,ok FROM ", db, " WHERE substance IN ('",
        slist,"') AND experimentator LIKE '",
        experimentator,"' AND ",testtype," LIKE '",
        type,"' AND ",
        whereClause," AND ok in (",
        ok,")",sep="")
    data <- sqlQuery(channel,query)
    odbcClose(channel)
    names(data)[[1]] <- "dose"
    names(data)[[2]] <- "response"
    data$dosefactor <- factor(data$dose)
    data$substance <- factor(data$substance,levels=substances)
    return(data)
}
    
drfit <- function(data, startlogEC50 = NA, chooseone=TRUE,
        lognorm = TRUE, logis = FALSE,
        linearlogis = FALSE, b0 = 2, f0 = 0)
{
    substances <- levels(data$substance)
    unit <- levels(as.factor(data$unit))
    logisf <- function(x,x0,b,k=1)
    {
        k / (1 + (x/x0)^b)
    }
    linearlogisf <- function(x,k,f,mu,b)
    {
        k*(1 + f*x) / (1 + ((2*f*(10^mu) + 1) * ((x/(10^mu))^b)))
    }

    ri <- rix <- 0                        # ri is the index over the result rows
                                          # rix is used later to check if any
                                          # model result was appended
    rsubstance <- array()                 # the substance names in the results
    rn <- vector()                        # number of dose-response curves 
    rlhd <- rlld <- vector()              # highest and lowest doses tested
    mtype <- array()                      # the modeltypes
    sigma <- array()                      # the standard deviation of the residuals
    logEC50 <- vector()
    stderrlogEC50 <- vector()
    slope <- vector()
    b <- vector()
    f <- vector()

    splitted <- split(data,data$substance)
    for (i in substances) {
        tmp <- splitted[[i]]
        fit <- FALSE
        n <- round(length(tmp$response)/9)
        if (length(tmp$response) == 0) {
            nodata = TRUE
        } else {
            nodata = FALSE
        }
        if (!nodata) {
            if (is.na(startlogEC50[i])){
                w <- 1/abs(tmp$response - 0.3)
                startlogEC50[[i]] <- sum(w * log10(tmp$dose))/sum(w)
            }
            highestdose <- max(tmp$dose)
            lowestdose <- min(tmp$dose)
            lhd <- log10(highestdose)
            lld <- log10(lowestdose)
            responseathighestdose <- mean(subset(tmp,dose==highestdose)$response)
            rix <- ri
            if (responseathighestdose < 0.5) {
                inactive <- FALSE

                if (linearlogis) {
                    m <- try(nls(response ~ linearlogisf(dose,1,f,logEC50,b),
                            data=tmp,
                            start=list(f=f0,logEC50=startlogEC50[[i]],b=b0)))
                    if (!inherits(m, "try-error")) {
                        fit <- TRUE
                        ri <- ri + 1
                        s <- summary(m)
                        sigma[[ri]] <- s$sigma
                        rsubstance[[ri]] <- i
                        rn[[ri]] <- n
                        rlld[[ri]] <- log10(lowestdose)
                        rlhd[[ri]] <- log10(highestdose)
                        mtype[[ri]] <- "linearlogis"
                        logEC50[[ri]] <- coef(m)[["logEC50"]]
                        slope[[ri]] <- NA
                        if (logEC50[[ri]] > rlhd[[ri]]) {
                            logEC50[[ri]] <- NA
                            stderrlogEC50[[ri]] <- NA
                            b[[ri]] <- NA
                            f[[ri]] <- NA
                        } else {
                            stderrlogEC50[[ri]] <- s$parameters["logEC50","Std. Error"]
                            b[[ri]] <- coef(m)[["b"]]
                            f[[ri]] <- coef(m)[["f"]]
                        }
                    }
                }

                if (logis) {
                # Instead of plogis(), the function logisf() defined above
                # could be used for regression against dose, not log10(dose)
                    m <- try(nls(response ~ plogis(-log10(dose),-logEC50,slope),
                            data=tmp,
                            start=list(logEC50=startlogEC50[[i]],slope=1)))
                    if (chooseone==FALSE || fit==FALSE) {
                        if (!inherits(m, "try-error")) {
                            fit <- TRUE
                            ri <- ri + 1
                            s <- summary(m)
                            sigma[[ri]] <- s$sigma
                            rsubstance[[ri]] <- i
                            rn[[ri]] <- n
                            rlld[[ri]] <- log10(lowestdose)
                            rlhd[[ri]] <- log10(highestdose)
                            mtype[[ri]] <- "logis"
                            logEC50[[ri]] <- coef(m)[["logEC50"]]
                            b[[ri]] <- NA
                            f[[ri]] <- NA
                            if (logEC50[[ri]] > rlhd[[ri]]) {
                                logEC50[[ri]] <- NA
                                slope[[ri]] <- NA
                                stderrlogEC50[[ri]] <- NA
                            } else {
                                slope[[ri]] <- coef(m)[["slope"]]
                                stderrlogEC50[[ri]] <- s$parameters["logEC50","Std. Error"]
                            }
                        }
                    }
                }

                if (lognorm) {
                    m <- try(nls(response ~ pnorm(-log10(dose),-logEC50,slope),
                                data=tmp,
                                start=list(logEC50=startlogEC50[[i]],slope=1)))
                    if (chooseone==FALSE || fit==FALSE) {
                        if (!inherits(m, "try-error")) {
                            fit <- TRUE
                            ri <- ri + 1
                            s <- summary(m)
                            sigma[[ri]] <- s$sigma
                            rsubstance[[ri]] <- i
                            rn[[ri]] <- n
                            rlld[[ri]] <- log10(lowestdose)
                            rlhd[[ri]] <- log10(highestdose)
                            mtype[[ri]] <- "lognorm"
                            logEC50[[ri]] <- coef(m)[["logEC50"]]
                            b[[ri]] <- NA
                            f[[ri]] <- NA
                            if (logEC50[[ri]] > rlhd[[ri]]) {
                                logEC50[[ri]] <- NA
                                slope[[ri]] <- NA
                                stderrlogEC50[[ri]] <- NA
                            } else {
                                slope[[ri]] <- coef(m)[["slope"]]
                                stderrlogEC50[[ri]] <- s$parameters["logEC50","Std. Error"]
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
            rn[[ri]] <- n
            if (nodata) {
                rlld[[ri]] <- rlhd[[i]] <- NA
                mtype[[ri]] <- "no data"
            } else {
                rlld[[ri]] <- log10(lowestdose)
                rlhd[[i]] <- log10(highestdose)
                if (inactive) {
                    mtype[[ri]] <- "inactive"
                } else {
                    mtype[[ri]] <- "no fit"
                }
            }
            sigma[[ri]] <- NA
            logEC50[[ri]] <- NA
            stderrlogEC50[[ri]] <- NA
            slope[[ri]] <- NA
            b[[ri]] <- NA
            f[[ri]] <- NA
        }
    }
    results <- data.frame(rsubstance,rn, rlhd, mtype, logEC50, stderrlogEC50, unit, sigma)
    names(results) <- c("Substance","n", "lhd","mtype","logEC50","std","unit","sigma")
    if (lognorm || logis) {
        results$slope <- slope
    }
    if (linearlogis) {
        results$b <- b
        results$f <- f
    }
    return(results)
}

drplot <- function(drresults, data = FALSE, dtype = "std", alpha = 0.95,
        path = "./", fileprefix = "drplot", overlay = FALSE,
        postscript = FALSE, 
        color = TRUE,
        datacolors = 1:8, fitcolors = "default")
{
    # Prepare plots
    devices <- 1
    if (postscript && !overlay) psdevices <- vector()
    if (!postscript && !overlay) x11devices <- vector()

    unit <- levels(as.factor(drresults$unit))

    # Get the plot limits on the x-axis (log of the dose)
    if(is.data.frame(data))
    {
        lld <- log10(min(data$dose))
        lhd <- log10(max(data$dose))
        hr <- max(data$response)
        dsubstances <- levels(data$substance)    
    } else {
        lld <- min(drresults[["logEC50"]],na.rm=TRUE) - 2
        lhd <- max(drresults[["logEC50"]],na.rm=TRUE) + 2
        if (length(subset(drresults,mtype=="linearlogis")$substance) != 0) {
            hr <- 1.8 
        } else {
            hr <- 1.0
        }
    }

    # Prepare overlay plot if requested
    if (overlay)
    {
        devices <- devices + 1
        if (postscript) {
            filename = paste(path,fileprefix,".eps",sep="")
            postscript(file=filename,
                    paper="special",width=7,height=7,horizontal=FALSE,pointsize=12) 
            cat("Created File: ",filename,"\n")
        } else {
            x11(width=7,height=7,pointsize=12)
        }
            
        plot(0,type="n",
            xlim=c(lld - 0.5, lhd + 2),
            ylim= c(-0.1, hr + 0.2),
            xlab=paste("Decadic Logarithm of the dose in ", unit),    
            ylab="Normalized response")
    }

    # Plot the data either as raw data or as error bars
    if(is.data.frame(data))
    {
        splitted <- split(data,data$substance)
        n <- 0
        for (i in dsubstances)
        {
            # Prepare the single graphs if an overlay is not requested
            if (!overlay)
            {
                devices <- devices + 1
                if (postscript) {
                    filename = paste(path,fileprefix,sub(" ","_",gsub("([\(\) ])", "", i)),".eps",sep="")
                    postscript(file=filename,
                            paper="special",width=7,height=7,horizontal=FALSE,pointsize=12)
                    psdevices[[i]] <- devices
                    cat("Created File: ",filename,"\n")
                } else {
                    x11(width=7,height=7,pointsize=12)
                    x11devices[[i]] <- devices
                }
                    
                plot(0,type="n",
                    xlim=c(lld - 0.5, lhd + 2),
                    ylim= c(-0.1, hr + 0.2),
                    xlab=paste("Decadic Logarithm of the dose in ", unit),    
                    ylab="Normalized response")
            }
            if (color == FALSE) datacolors <- rep("black",length(dsubstances))
            n <- n + 1
            if (!overlay) legend(lhd - 1, hr + 0.1, i,lty = 1, col = datacolors[[n]])
            tmp <- splitted[[i]]
            tmp$dosefactor <- factor(tmp$dose)  # necessary because the old
                                                # factor has all levels, not 
                                                # only the ones tested with
                                                # this substance
            if (dtype == "raw") {
                points(log10(tmp$dose),tmp$response,col=datacolors[[n]])
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
                segments(x,bottoms,x,tops,col=datacolors[[n]])
                points(x,means,col=datacolors[[n]])
                smidge <- 0.05
                segments(x - smidge,bottoms,x + smidge,bottoms,col=datacolors[[n]])
                segments(x - smidge,tops,x + smidge,tops,col=datacolors[[n]])
            }
        }
    }

    # Plot the fitted dose response models from drresults
    fits <- subset(drresults,!is.na(logEC50))
    nf <-  length(fits$Substance)  # number of fits to plot
    if (fitcolors[[1]] == "default")
    {
        if (nf <= 8)
        {
            defaultfitcolors <- palette()
        } else
        {
            defaultfitcolors <- rainbow(nf)
        }
    }
    legendcolors <- vector()
    for (i in 1:nf)
    {
        s <- as.character(fits[i,"Substance"]) # The substance name to display
        if (!overlay && !is.data.frame(data))
        {
            devices <- devices + 1
            if (postscript) {
                filename = paste(path,fileprefix,sub(" ","_",gsub("([\(\) ])", "", s)),".eps",sep="")
                postscript(file=filename,
                        paper="special",width=7,height=7,horizontal=FALSE,pointsize=12) 
                psdevices[[s]] <- devices
                cat("Created File: ",filename,"\n")
            } else {
                x11(width=7,height=7,pointsize=12)
                x11devices[[s]] <- devices
            }
                
            plot(0,type="n",
                xlim=c(lld - 0.5, lhd + 2),
                ylim= c(-0.1, hr + 0.2),
                xlab=paste("Decadic Logarithm of the dose in ", unit),    
                ylab="Normalized response")
        }
        if (postscript && !overlay) {
            dev.set(psdevices[[s]]) }
        if (!postscript && !overlay) {
            dev.set(x11devices[[s]]) }

        if (color == FALSE) {
            fitcolor <- "black" 
        } else {
            if (fitcolors[[1]] == "default")
            {
                fitcolor <- defaultfitcolors[[i]]
            } else {
                fitcolor <- fitcolors[[i]]
            }
        }
        if (!overlay) legend(lhd - 1, hr + 0.1, s,lty = 1, col = fitcolor)
        legendcolors[[i]] <- fitcolor
        logEC50 <- fits[i,"logEC50"]
        mtype <- as.character(fits[i, "mtype"])
        if (mtype == "lognorm")
        {
            slope <- fits[i,"slope"]
            plot(function(x) pnorm(-x,-logEC50,slope),lld - 0.5, lhd + 2, add=TRUE,col=fitcolor)
        }
        if (mtype == "logis")
        {
            slope <- fits[i,"slope"]
            plot(function(x) plogis(-x,-logEC50,slope),lld - 0.5, lhd + 2, add=TRUE,col=fitcolor)
        }
    }
    if (overlay) {
        legend(lhd - 1, hr + 0.1, as.vector(fits$Substance), lty = 1, col = legendcolors)
    }
    if (devices > 1 && postscript)
    {
        for (i in 2:devices) {
            dev.off(i)
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
        
        cat("Plate ",plate," from database ",db,"\n",
            "\tExperimentator: ",levels(platedata$experimentator),"\n",
            "\tType(s): ",levels(platedata$type),"\n",
            "\tPerformed on : ",levels(platedata$performed),"\n",
            "\tSubstance(s): ",levels(platedata$substance),"\n",
            "\tConcentration unit: ",levels(platedata$unit),"\n",
            "\tOK: ",levels(platedata$ok),"\n",
            "\t\tNumber \tMean \tStandard Deviation\n",
            "blind\t\t",numberOfBlinds,"\t",meanOfBlinds,"\t",stdOfBlinds,"\n",
            "control\t",numberOfControls,"\t",meanOfControls,"\t",stdOfControls,"\n")
        
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

        legend(3.0,1.5,substances, pch=1, col=1:length(substances))
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
    
    legend(3.5,1.7,plates, pch=1, col=1:length(plates))
    title(main=paste(substance," - ",levels(data$experimentator)," - ",levels(data$type)))
 
    cat("Substanz ",substance,"\n",
        "\tExperimentator(s):",levels(data$experimentator),"\n",
        "\tType(s):\t",levels(data$type),"\n",
        "\tSubstance(s):\t",levels(data$substance),"\n",
        "\tPlate(s):\t",plates,"\n\n")
}
