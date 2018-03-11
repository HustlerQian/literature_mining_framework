
xmlTovec2 <- function(pubmed_data) {
  options(warn = -1)
  tmp.article <- custom_grep(xml_data = pubmed_data, tag = "PubmedArticle", format = "char")
  
  tmp.abstract <- NA
  tmp.date <- c(NA,NA,NA)
  tmp.jabbrv <- NA
  tmp.title <- custom_grep(xml_data = tmp.article, tag = "ArticleTitle", format = "char")
  if(is_empty(tmp.title))tmp.title <- NA
  tmp.abstract <- custom_grep(xml_data = tmp.article, tag = "AbstractText", format = "char")
  if (length(tmp.abstract) > 1){
    tmp.abstract <- paste(tmp.abstract, collapse = " ", sep = " ")
  } else if (length(tmp.abstract) < 1) {
    tmp.abstract <- NA
  }
  tmp.date <- custom_grep(xml_data = tmp.article, tag = "PubDate", format = "char")
  tmp.date <- sapply(c("Year", "Month", "Day"), (function(tt){
    custom_grep(xml_data = tmp.date, tag = tt, format = "char")   
  }))
  tmp.jabbrv  <- custom_grep(xml_data = tmp.article, tag = "ISOAbbreviation", format = "char")
  
  tmp.resout <- c(title=tmp.title,
                  abstract=tmp.abstract,
                  year = as.character(tmp.date[1]),
                  month = as.character(tmp.date[2]),
                  day = as.character(tmp.date[3]),
                  jabbrv=tmp.jabbrv)
  options(warn = 0)
  return(tmp.resout)
}
