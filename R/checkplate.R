checkplate <- function(id, db = c("cytotox", "enzymes"))
{
  db <- match.arg(db)
  checkexperiment(id, db = db)
}
