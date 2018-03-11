getabs <- function(list) {
  if (length(list) > 2) {
    abstract <- paste(list[1:length(list) - 1], collapse = " ", 
                      sep = " ")
  }
  else if (length(list) < 1) {
    abstract <- 'NA'
  }
  else {
    abstract <- list[1]
  }
  return(abstract)
}
record <- xml2 %>% xml2::xml_find_all(".//PubmedArticle//MedlineCitation//Article")
journal <- record %>% xml2::xml_find_all("Journal//ISOAbbreviation") %>% 
  xml2::xml_contents() %>% as.character() %>% unlist()
title <- record %>% xml2::xml_find_all("ArticleTitle") %>% 
  xml2::xml_contents() %>% as.character() %>% unlist()
year <- record %>% xml2::xml_find_all("Journal//JournalIssue//PubDate") %>% 
  xml2::as_list() %>% purrr::map(unlist) %>% purrr::map(paste0) %>% 
  purrr::map(`[`, 1) %>% unlist()
month <- record %>% xml2::xml_find_all("Journal//JournalIssue//PubDate") %>% 
  xml2::as_list() %>% purrr::map(unlist) %>% purrr::map(paste0) %>% 
  purrr::map(`[`, 2) %>% unlist()
day <- record %>% xml2::xml_find_all("Journal//JournalIssue//PubDate") %>% 
  xml2::as_list() %>% purrr::map(unlist) %>% purrr::map(paste0) %>% 
  purrr::map(`[`, 3) %>% unlist()



abstract <- record %>% xml2::as_list() %>% purrr::map(`[[`, 
  "Abstract") %>% purrr::map(getabs) %>% purrr::map_chr(unlist) %>% 
  stringr::str_replace_all("list\\(\"", "") %>%
  stringr::str_replace_all("\"\\)", "") %>%
  unlist() 

paperdf <- tibble::as.tibble(cbind(journal, title, year,month, day, abstract)) %>% 
  dplyr::mutate(day = dplyr::case_when(is.na(day) ~"01", !is.na(day) ~ day)) %>% 
  dplyr::mutate(month = dplyr::case_when(!(month %in% month.abb) & is.na(month) ~ "Jan", !(month %in% month.abb) & !is.na(month) ~ month.abb[as.numeric(month)], month %in% month.abb ~ month)) %>%
  tidyr::unite(date, year, month, day, sep = "") %>%
  dplyr::mutate(date = as.Date(date, "%Y%b%d")) %>%
  dplyr::bind_cols(line = 1:length(title))
                                                                                                     



options(warn = 0)
return(paperdf)