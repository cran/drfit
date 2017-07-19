if(getRversion() >= '2.15.1') utils::globalVariables(c("type", "conc", "substance"))
checkexperiment <- function(id, db = "ecotox", endpoint = "%")
{
    databases <- data.frame(
        responsename=c("viability","activity","raw_response"),
        testtype=c("celltype","enzyme","organism"),
        exptype=c("plate","plate","experiment"))
    rownames(databases) <- c("cytotox","enzymes","ecotox")

    if (!(db %in% rownames(databases))) stop("Database is not supported")

    if (requireNamespace("RODBC")) {
      channel <- RODBC::odbcConnect(db, uid="cytotox", pwd="cytotox", case="tolower")
    } else {
      stop("For this function, the RODBC package has to be installed and configured.")
    }

    responsename = as.character(databases[db,1])
    testtype = as.character(databases[db,2])
    exptype = as.character(databases[db,3])

    exptable <- paste(exptype, "s", sep="")
    commentquery <- paste("SELECT comment FROM ", exptable ,
        " WHERE ", exptype, " = ", id)
    commentdata <- RODBC::sqlQuery(channel,commentquery)
    comment <- as.character(commentdata[[1]])

    expquery <- paste("SELECT experimentator, substance, ",
        testtype, ", conc, unit,", responsename, ", type, raw_0, duration, performed, ok",
        " FROM ",db," WHERE ",exptype,"=", id,
            sep = "")

    if (db == "ecotox") {
        expquery <- paste0(expquery, " AND type LIKE '", endpoint, "'") 
    }

    expdata <- RODBC::sqlQuery(channel,expquery)

    if (db %in% c("cytotox","enzymes")) {
        controlquery <- paste("SELECT type,response FROM controls
            WHERE plate=",id)
        controldata <- RODBC::sqlQuery(channel,controlquery)
    }

    RODBC::odbcClose(channel)

    op <- par(ask=TRUE)
    on.exit(par(op))

    if (db %in% c("cytotox","enzymes")) {
        blinds <- subset(controldata, type == "blind")
        controls <- subset(controldata, type == "control")
        QA <- matrix(nrow = 2, ncol = 4,
            dimnames = list(c("Blind", "Control (conc = 0)"),
                            c("Number", "Mean", "Std. Dev.", "% Std. Dev")))

        QA[1, 1] <- length(blinds$response)
        QA[1, 2] <- signif(mean(blinds$response), 2)
        QA[1, 3] <- signif(sd(blinds$response), 2)
        QA[1, 4] <-signif(QA[1, 3] * 100 / QA[1, 2],2)
    } else {
        # Use raw response for ecotox
        expdata$response <- expdata$raw_response
        
        if (nlevels(expdata$type) > 1) {
          message("There are data for more than one type of raw response in your data.\n",
                  "The types are ", paste(levels(expdata$type), collapse = " and "), ".\n",
                  "You should choose one of these types using 'endpoint = \"$type\"'",
                  "in your call to checkexperiment\n",
                  "For now, we are continuing with the data for ", levels(expdata$type)[1])
        }
        endpoint <- expdata$type[1]
        expdata <- subset(expdata, type == endpoint)

        controls <- subset(expdata, conc == 0)
        expdata <- subset(expdata, conc != 0)

        QA <- matrix(nrow = 1, ncol = 4,
            dimnames = list(c("Control (conc = 0)"),
                            c("Number", "Mean", "Std. Dev.", "% Std. Dev")))

    }

    numberOfControls <- length(controls$response)
    QA["Control (conc = 0)", 1] <- numberOfControls
    if (numberOfControls > 0) {
        QA["Control (conc = 0)", 2] <- signif(mean(controls$response),2)
        QA["Control (conc = 0)", 3] <- signif(sd(controls$response),2)
        QA["Control (conc = 0)", 4] <- signif(QA["Control (conc = 0)", 3] * 100 /
                                              QA["Control (conc = 0)", 2],2)
    }

    if (db == "ecotox") {
        if (identical(as.character(levels(expdata$organism)), "Vibrio fischeri")) {
            positive <- subset(expdata, substance == "Na Cl")
            if (nrow(positive) > 0) {
                 QA <- rbind(QA,
                             c(nrow(positive),
                               signif(mean(positive$raw_response), 2),
                               signif(sd(positive$raw_response), 2),
                               signif(100 * sd(positive$raw_response) /
                                      mean(positive$raw_response), 2)))

                 rownames(QA) <- c("Control (conc = 0)",
                                   "Positive control (Na Cl)")
            }
            expdata <- subset(expdata, substance != "Na Cl", drop = TRUE)
        }
    }

    if (length(expdata$experimentator) < 1) {
        stop("There is no response data for ",exptype," ",
            id," in database ",db,"\n")
    }
    exptypestring <- paste(toupper(substring(exptype,1,1)),
        substring(exptype,2),sep="")
    expdata$experimentator <- factor(expdata$experimentator)
    expdata$type <- factor(expdata[[testtype]])
    expdata$performed <- factor(as.character(expdata$performed))
    expdata$substance <- factor(expdata$substance)
    expdata$unit <- factor(expdata$unit)
    expdata$ok <- factor(expdata$ok)

    # Info on the experiment
    cat("\n",exptypestring,id,"from database",db,":\n\n",
        "\tExperimentator(s):\t",levels(expdata$experimentator),"\n",
        "\tType(s):\t\t",levels(expdata$type),"\n",
        "\tPerformed on:\t\t",levels(expdata$performed),"\n",
        "\tSubstance(s):\t\t",levels(expdata$substance),"\n",
        "\tConcentration unit(s):\t",levels(expdata$unit),"\n",
        "\tComment:\t\t",comment,"\n",
        "\tOK Levels:\t\t",levels(expdata$ok),"\n\n")

    print(QA)

    # Control growth rate for Lemna and algae
    if (endpoint %in% c("cell count", "frond area", "frond number")) {
      duration <- unique(expdata$duration) # in hours
      if (length(duration) > 1) stop("More than one duration in the data")
      response_0 <- unique(expdata$raw_0)
      if (length(response_0) > 1) stop("More than one mean response at time 0 in the data")
      t_days <- duration / 24
      control_growth_rates <- (log(controls$response) - log(response_0)) / t_days
      cat("\nMean growth rate in controls:\t", round(mean(control_growth_rates), 3), "per day\n")
    }


    # Box plot of control data
    if (db == "ecotox") {
        boxplot(controls$response,
            names="controls",
            ylab=endpoint,
            ylim=range(controls$response, na.rm = TRUE),
            boxwex=0.4,
            main=paste("Plate ",id))
    } else {
        boxplot(blinds$response,controls$response,
            names=c("blinds","controls"),
            ylab="Response",
            boxwex=0.4,
            main=paste("Plate ",id))
    }

    # Plot of dose response data
    drdata <- expdata[c(2,4,6)]
    drdata$substance <- factor(drdata$substance)
    substances <- levels(drdata$substance)

    lld <- log10(min(subset(drdata,conc!=0)$conc))
    lhd <- log10(max(drdata$conc))

    ylab <- if (db == "ecotox") endpoint
            else responsename

    plot(1,type="n",
        xlim = c(lld - 0.5, lhd + 2),
        ylim = range(expdata[responsename], na.rm = TRUE),
        xlab = paste("decadic logarithm of the concentration in ",levels(expdata$unit)),
        ylab = ylab)

    drdatalist <- split(drdata,drdata$substance)

    for (i in 1:length(drdatalist)) {
        points(log10(drdatalist[[i]]$conc),drdatalist[[i]][[responsename]],col=i);
    }

    legend("topright",substances, pch=1, col=1:length(substances), inset=0.05)
    title(main=paste(levels(expdata$experimentator),
        " - ",levels(expdata$type)))
}
