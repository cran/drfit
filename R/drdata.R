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
    data$substance <- factor(data$substance,levels=substances)
    return(data)
}
