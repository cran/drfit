checkexperiment <- function(id,db="ecotox")
{
    databases <- data.frame(
        responsetype=c("viability","activity","response"),
        testtype=c("celltype","enzyme","organism"),
        exptype=c("plate","plate","experiment"))
    rownames(databases) <- c("cytotox","enzymes","ecotox")

    if (!(db %in% rownames(databases))) stop("Database does not exist")

    library(RODBC) 
    channel <- odbcConnect(db,uid="cytotox",pwd="cytotox",case="tolower")

    responsetype = as.character(databases[db,1])
    testtype = as.character(databases[db,2])
    exptype = as.character(databases[db,3])
        
    expquery <- paste("SELECT experimentator,substance,",
        testtype,",conc,unit,",responsetype,",performed,ok",
        "FROM ",db," WHERE ",exptype,"=", id)

    expdata <- sqlQuery(channel,expquery)

    if (db %in% c("cytotox","enzymes")) {
        controlquery <- paste("SELECT type,response FROM controls 
            WHERE plate=",id)
        controldata <- sqlQuery(channel,controlquery)
    }

    odbcClose(channel)

    op <- par(ask=TRUE)
    on.exit(par(op))

    if (db %in% c("cytotox","enzymes")) {
        blinds <- subset(controldata,type=="blind")
        controls <- subset(controldata,type=="control")

        numberOfBlinds <- length(blinds$response)
        meanOfBlinds <- signif(mean(blinds$response),2)
        stdOfBlinds <- signif(sd(blinds$response),2)
    } else {
        controls <- subset(expdata,conc == 0)
        expdata <- subset(expdata, conc != 0)

        numberOfBlinds <- NA
        meanOfBlinds <- NA
        stdOfBlinds <- NA
        
    }
    numberOfControls <- length(controls$response)
    if (numberOfControls > 0) {
        meanOfControls <- signif(mean(controls$response),2)
        stdOfControls <- signif(sd(controls$response),2)
        percentstdOfcontrols <-signif(stdOfControls *100/meanOfControls,2)
    } else {
        meanOfControls <- stdOfControls <- percentstdOfcontrols <- NA
    }

    if (length(expdata$experimentator) < 1) {
        stop("There is no response data for ",exptype," ",
            id," in database ",db,"\n")
    } 
    expdata$experimentator <- factor(expdata$experimentator)
    expdata$type <- factor(expdata[[testtype]])
    expdata$performed <- factor(as.character(expdata$performed))
    expdata$substance <- factor(expdata$substance)
    expdata$unit <- factor(expdata$unit)
    expdata$ok <- factor(expdata$ok)
    
    cat("\n",exptype," ",id," from database ",db,"\n",
        "\tExperimentator(s): ",levels(expdata$experimentator),"\n",
        "\tType(s): ",levels(expdata$type),"\n",
        "\tPerformed on: ",levels(expdata$performed),"\n",
        "\tSubstance(s): ",levels(expdata$substance),"\n",
        "\tConcentration unit(s): ",levels(expdata$unit),"\n",
        "\tOK: ",levels(expdata$ok),"\n",
        "\t\tNumber \tMean \tStd. Dev. \t% Std. Dev.\n",
        "\tblind\t",numberOfBlinds,"\t",meanOfBlinds,"\t",stdOfBlinds,"\n",
        "\tcontrol\t",numberOfControls,"\t",meanOfControls,"\t",
            stdOfControls,"\t\t",percentstdOfcontrols,"\n")
    

    if (db == "ecotox") {
        boxplot(controls$response,
            names="controls",
            ylab="Response",
            ylim=c(0,max(controls$response)),
            boxwex=0.4,
            main=paste("Plate ",id))
    } else {
        boxplot(blinds$response,controls$response,
            names=c("blinds","controls"),
            ylab="Response",
            boxwex=0.4,
            main=paste("Plate ",id))
    }
    
    drdata <- expdata[c(2,4,6)]
    drdata$substance <- factor(drdata$substance)
    substances <- levels(drdata$substance)
       
    lld <- log10(min(subset(drdata,conc!=0)$conc))
    lhd <- log10(max(drdata$conc))

    plot(1,type="n",
        xlim=c(lld - 0.5, lhd + 2), 
        ylim= c(-0.1, 2), 
        xlab=paste("decadic logarithm of the concentration in ",levels(expdata$unit)),
        ylab=responsetype)
    
    drdatalist <- split(drdata,drdata$substance)
    
    for (i in 1:length(drdatalist)) {
        points(log10(drdatalist[[i]]$conc),drdatalist[[i]][[responsetype]],col=i);
    }

    legend("topright",substances, pch=1, col=1:length(substances), inset=0.05)
    title(main=paste(levels(expdata$experimentator),
        " - ",levels(expdata$type)))
}
