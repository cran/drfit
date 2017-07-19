if(getRversion() >= '2.15.1') utils::globalVariables(c("type", "conc"))
checkcontrols <- function(last = 10, id = NULL, db = "cytotox",
                          celltype = "IPC-81", enzymetype = "AChE",
                          organism = "Vibrio fischeri",
                          endpoint = "%", qcc = c("R", "xbar"))
{

    if (!(db %in% c("cytotox", "ecotox", "enzymes"))) stop("Database is not supported")

    if (requireNamespace("RODBC")) {
      channel <- RODBC::odbcConnect(db, uid="cytotox", pwd="cytotox", case="tolower")
    } else {
      stop("For this function, the RODBC package has to be installed and configured.")
    }

    if (db %in% c("cytotox","enzymes")) {
        if (is.null(id[1])) {
            platequery <- "SELECT plate FROM"
            if (db == "cytotox") {
                platequery <- paste0(platequery, " cytotox WHERE celltype like
                                     '", celltype, "'")
            }
            if (db == "enzymes") {
                platequery <- paste0(platequery, " enzymes WHERE enzyme like
                                     '", enzymetype, "'")
            }
                platequery <- paste(platequery,
                    "GROUP BY plate ORDER BY plate DESC LIMIT", last)
            plates <- RODBC::sqlQuery(channel, platequery)$plate
        } else {
            plates <- id
        }
        controlquery <- paste("SELECT plate, type, location, response FROM controls",
            "WHERE plate IN (", paste(plates, collapse = ", "), ")")
        controldata <- RODBC::sqlQuery(channel,controlquery)
    } else {
        if (is.null(id[1])) {
            lastquery = paste0("SELECT experiment, type FROM ecotox ",
                "WHERE organism LIKE '", organism, "'",
                "AND type LIKE '", endpoint, "' ",
                "GROUP BY experiment ORDER BY experiment DESC LIMIT ", last)
            res <- RODBC::sqlQuery(channel, lastquery)
            if (nrow(res) == 0) {
                stop("No results for endpoint", endpoint)
            } else {
                if (nlevels(res$type) > 1) {
                    stop("Found more than one endpoint type:\n",
                         paste(levels(res$type), collapse = ", "), "\n",
                         "Please specify an endpoint in your call to checkcontrols()")
                }
                experiments <- res$experiment
            }
        } else {
            experiments <- id
        }
        expquery <- paste0("SELECT ",
            "experimentator, experiment, substance, organism, type, conc, unit, raw_response, ",
            "performed, ok ",
            "FROM ecotox ",
            "WHERE experiment IN (", paste(experiments, collapse = ", "), ") ",
            "AND organism LIKE '", organism, "' ",
            "AND type LIKE '", endpoint, "'")
        expdata <- RODBC::sqlQuery(channel, expquery)
        if (nlevels(expdata$type) > 1) {
            stop("Found more than one endpoint type:\n",
                 paste(levels(expdata$type), collapse = ", "), "\n",
                 "Please specify an endpoint in your call to checkcontrols()")
        }
        # Use the raw response for QA
        expdata$response <- expdata$raw_response
    }

    RODBC::odbcClose(channel)

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
        controls <- subset(expdata, conc == 0)
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

    # The report
    cat("\nDatabase", db, "\n")

    if (db == "ecotox") {
        cat("Organism", organism, "\n")
        cat("Endpoint", unique(expdata$type), "\n")
        cat("\nExperiments ", paste(experiments, collapse = ", "), "\n\n")
    } else {
        if (db == "cytotox") cat ("Cell type", celltype, "\n")
        if (db == "enzymes") cat ("Enzyme type", enzymetype, "\n")
        cat("\nPlates", paste(plates, collapse = ", "), "\n\n")
    }
    print(QA)

    if (!is.na(qcc[1])) {
        op <- par(ask=TRUE)
        on.exit(par(op))
        requireNamespace("reshape2")
        if (db == "ecotox") {
            controls$row <- rownames(controls)
            controls_molten <- melt(controls[c("experiment", "row", "response")],
                                    id = c("experiment", "row"))
            controls_wide <- acast(controls_molten, formula = experiment ~ row)

        } else {
            controls_molten <- melt(controls[c("plate", "location", "response")],
                                    id = c("plate", "location"))
            controls_wide <- acast(controls_molten, formula = plate ~ location)
        }
        if ("R" %in% qcc) {
            qcc(controls_wide, type = "R", nsigmas = 3,
                title = "Range chart",
                data.name = "Controls (conc = 0)")
        }
        if ("xbar" %in% qcc) {
            qcc(controls_wide, type = "xbar", nsigmas = 3,
                title = "Mean chart",
                data.name = "Controls (conc = 0)")
        }
    }

    invisible(controls)
}
# vim: ts=4 sw=4 expandtab:
