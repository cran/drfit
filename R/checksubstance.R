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
