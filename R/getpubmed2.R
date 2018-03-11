getpubmed2 <-function (query, start = 1, end = 100) 
{
query <- as.character(query)
query <- gsub(" ", "+", query, fixed = TRUE)
PID <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?", 
             "db=pubmed&term=", query, "&usehistory=y", sep = "")
xml <- xml2::read_xml(PID)
list <- xml2::as_list(xml)$eSearchResult
n <- list$Count[[1]]
warning <- paste(n, "records founds")
message(warning)
efetch_url = paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?", 
                   "db=pubmed&WebEnv=", list$WebEnv[[1]], "&query_key=", 
                   list$QueryKey[[1]], "&retstart=", start - 1, "&retmax=", 
                   end, "&retmode=xml", sep = "")


xml2 <- xml2::read_xml(efetch_url)
return(xml2)
}
