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

        boxplot(blinds$response,controls$response,
            names=c("blinds","controls"),
            ylab="Response",main=paste("Plate ",plate))
        
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
