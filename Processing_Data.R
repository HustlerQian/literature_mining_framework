library(easyPubMed)
library(tidyverse)
#library(xml2)
#library(lubridate)
#library(scifetch)
rm(list = ls())
#Use functions 
source('./R/getpubmed2.R')
source('./R/xmlTovec2.R')



path_data <- './cas9'
#batch_pubmed_download
if(!dir.exists(path_data)){dir.create(path_data)}
query <- 'crispr[All]'
# '"clustered regularly interspaced short palindromic repeats"[MeSH Terms] OR 
# ("clustered"[All Fields] AND "regularly"[All Fields] 
# AND "interspaced"[All Fields] AND "short"[All Fields] 
# AND "palindromic"[All Fields] AND "repeats"[All Fields]) OR 
# "clustered regularly interspaced short palindromic repeats"[All Fields] OR
# "crispr"[All Fields]'

#query <- 'Solid Phase MicroExtraction[MH] AND 1997/08:2017/08[DP]'

#Downloading
# out <- batch_pubmed_download(pubmed_query_string = query, 
#                              dest_file_prefix = "raw", batch_size = 200, 
#                              dest_dir = './cas9/')
out <- list.files(path = './cas9/', pattern = "*.xml",full.names = T)

paperdf <- tibble()
for (i in 1:length(out)){
  xmlfile <- xmlParse(out[i])
  paper.data <- articles_to_list(xmlfile)
  papers.list <- t(sapply(1:length(paper.data), (function(i){
    out <-  tryCatch(xmlTovec2(paper.data[i]), error = function(e) {NULL})
  })))
  #Use batch size for debug
  papertemp <- tibble()
  if(dim(papers.list)[1]==200){
  papertemp <- as_tibble(papers.list)
  }
  else{
    for( i in 1:length(papers.list)){
      temp <- data.frame(
        title = papers.list[[i]]['title'],
        abstract = papers.list[[i]]['abstract'],
        year = papers.list[[i]]['year'],
        month = papers.list[[i]]['month'],
        day = papers.list[[i]]['day'],
        jabbrv = papers.list[[i]]['jabbrv']
      )
    }
    papertemp <- bind_rows(papertemp,temp)
  }
  paperdf <- bind_rows(paperdf,papertemp)
}
#Exclude NA
paperdf_exna<-paperdf[!is.na(paperdf$title),]
paperdf_exna<-paperdf_exna[!is.na(paperdf_exna$abstract),]
paperdf_exna<-paperdf_exna[!is.na(paperdf_exna$year),]
paperdf_exna<-paperdf_exna[!is.na(paperdf_exna$jabbrv),]

tmdf <- paperdf_exna
#?save
save(tmdf,file = 'cas9.Rdata')
#tmdf <- getpubmed2(query, start = 1, end = 10000) 
  #%>%
  #getpubmedtbl() %>%
  #mutate(time = as.POSIXct(date, origin = "1970-01-01"),
  #       month = round_date(date, "month"))
