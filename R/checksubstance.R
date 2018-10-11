checksubstance <- function(substance, 
  db = c("cytotox", "enzymes", "ecotox"), 
  experimentator = "%",
  celltype = "%", enzymetype = "%", organism = "%",
  endpoint = "%",
  whereClause = "1", ok= "%")
{
  db = match.arg(db)

  databases <- data.frame(
    responsename = c("viability", "activity", "raw_response"),
    testtype = c("celltype", "enzyme", "organism"),
    exptype = c("plate", "plate", "experiment"),
    row.names = c("cytotox", "enzymes", "ecotox"),
    stringsAsFactors = FALSE)

  con <- dbConnect(odbc(), "cytotox", database = db)

  responsename <- databases[db, 1]
  testtype <- databases[db, 2]
  exptype <- databases[db, 3]

  if (db == "cytotox") {
      type <- celltype
  }
  if (db == "enzymes") {
      type <- enzymetype
  }
  if (db == "ecotox") {
      type <- organism
  }

  query <- paste("SELECT experimentator,substance,",
      testtype, ",", exptype, ",conc,unit,",responsename,",ok",
      " FROM ",db," WHERE substance LIKE '",
      substance,"' AND experimentator LIKE '",
      experimentator,"' AND ",testtype," LIKE '",
      type,"' AND ",
      whereClause," AND ok LIKE '",ok,"'",
              sep = "")

  if (db == "ecotox") {
      query <- paste(query, " AND type LIKE '",
              endpoint, "'", sep = "")
  }

  data <- dbGetQuery(con, query)

  if (length(data$experimentator) < 1) {
      stop(paste("\nNo response data for",substance,"in database",
              db,"found with these parameters\n"))
  }

  data$experimentator <- factor(data$experimentator)
  data$substance <- factor(data$substance)
  substances <- levels(data$substance)
  data$type <- factor(data[[testtype]])
  data[[exptype]] <- factor(data[[exptype]])
  experiments <- levels(data[[exptype]])
  concentrations <- split(data$conc,data$conc)
  concentrations <- as.numeric(names(concentrations))
  data$unit <- factor(data$unit)
  data$ok <- factor(data$ok)

  if (length(experiments)>6) {
      palette(rainbow(length(experiments)))
  }

  plot(log10(data$conc),data[[responsename]],
      xlim=c(-2.5, 4.5),
      ylim= range(data[[responsename]], na.rm = TRUE),
      xlab=paste("decadic logarithm of the concentration in ",levels(data$unit)),
      ylab=responsename)

  explist <- split(data,data[[exptype]])

  for (i in 1:length(explist)) {
      points(log10(explist[[i]]$conc),explist[[i]][[responsename]],col=i);
  }

  legend("topleft", experiments, pch=1, col=1:length(experiments), inset=0.05)
  title(main=paste(substance," - ",levels(data$experimentator)," - ",levels(data$type)))

  exptypename <-  paste(toupper(substring(exptype,1,1)),
      substring(exptype,2), sep = "")
  experimentators <- paste(levels(data$experimentator), collapse = " ")
  types <- paste(levels(data$type), collapse = " ")
  experiments <- paste(levels(data[[exptype]]), collapse = " ")
  class(experiments)
  cat("\n\tSubstanz:\t\t",substance,"\n",
      "\tExperimentator(s):\t", experimentators,"\n",
      "\tType(s): \t\t",types,"\n",
      "\tEndpoint: \t\t",endpoint,"\n",
      "\t", exptypename, "(s):\t\t",experiments,"\n\n", sep = "")
}
