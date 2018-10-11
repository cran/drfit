drdata <- function(substances, experimentator = "%", db = "cytotox",
    celltype = "IPC-81", enzymetype = "AChE",
    organism = "Vibrio fischeri", endpoint = "Luminescence", whereClause = "1",
    ok = "'ok','no fit'")
{
  # Connect to the correct database via the DSN
  con <- dbConnect(odbc(), "cytotox", database = db)

  # Construct the query
  slist <- paste(substances, collapse = "','")
  if (db == "cytotox") {
    experimenttype <- "plate"
    responsetype <- "viability"
    testtype <- "celltype"
    type <- celltype
  } else {
    if (db == "enzymes") {
      experimenttype <- "plate"
      responsetype <- "activity"
      testtype <- "enzyme"
      type <- enzymetype
    } else {
      experimenttype <- "experiment"
      responsetype <- "response"
      testtype <- "organism"
      type <- organism
    }
  }

  query <- paste0(
    "SELECT conc,", responsetype, ", unit, experimentator, ",
            experimenttype, ", substance, ", testtype,
            ", ok ",
      "FROM ", db, " ",
      "WHERE ",
        "substance IN ('", slist, "') AND ",
        "experimentator LIKE '", experimentator,"' AND ",
        testtype, " LIKE '", type, "' AND ",
        whereClause, " AND ",
        "ok in (", ok, ")")

  if (db == "ecotox") query <- paste0(query, " AND type LIKE '", endpoint, "'")

  # Get the data, format and return them
  data <- dbGetQuery(con, query)

  names(data)[[1]] <- "dose"
  names(data)[[2]] <- "response"
  data$substance <- factor(data$substance, levels = substances)

  return(data)
}
