#Library Section ---------------------------

library(tidyverse)
library(ggplot2)
library(tidytext)
library(stringr)
library(scales)
library(lubridate)
library(topicmodels)

options(stringsAsFactors = F)
cas9.df <- read.delim('D:/wyq/paperNetwork/data/result/rdata/conceptsInFile2_all.txt',sep = '\t')

colnames(cas9.df) <- c('PMID','Date','Concept','Count','SMtype','Category')

cas9.tdf <- tibble::as.tibble(cas9.df)

#Make 100 PMID for annotation
PMID100 <- sample(unique(cas9.tdf$PMID),100)
#Output PMID100
write(PMID100,file = 'annotation100.txt')
#Count the concepts in different category
cas9.tdf_100 <- cas9.tdf %>%
  filter(PMID %in% PMID100)
table(cas9.tdf_100$Category)

#PMID 
#Filter with gene,species and disease
cas9.tdf <- cas9.tdf %>%
  filter(Category %in% c('Disease_pathologic','Material_gene','Material_species'))

disease_name <- cas9.tdf[cas9.tdf$Category=='Disease_pathologic',]$Concept
gene_name <- cas9.tdf[cas9.tdf$Category=='Material_gene',]$Concept
species_name <- cas9.tdf[cas9.tdf$Category=='Material_species',]$Concept

#Filter by concept using top10 tfidf score
tf_idf <- cas9.tdf %>%
  bind_tf_idf(PMID, Concept, Count) %>%
  arrange(desc(tf_idf))

tf_idf_top10 <- tf_idf %>%
  group_by(PMID) %>%
  #filter(PMID==25775608)
  arrange(desc(PMID))%>%
  top_n(10)

#To see gene and disease relationship by PMI
abs_word_pairs_PMI <- tf_idf_top10 %>%
  #filter(Category %in% c('Material_gene','Material_species')) %>%
  pairwise_pmi(Concept,PMID,sort = TRUE) %>%
  filter(item1 %in% species_name) %>%
  #filter(item1 %in% gene_name) %>%
  filter(item2 %in% disease_name)


abs_word_pairs_PMI %>%
  head(100) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = pmi, edge_width = pmi), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "PMI Concepts in abstract") +
  theme_void()


abs_word_pairs_cor <- cas9.tdf %>%
  #filter(Concept %in% index) %>%
  #filter(Count >5) %>%
  pairwise_cor(Concept,PMID,Count)%>%
  arrange(desc(correlation)) 

abs_word_pairs_cor %>%
  filter(item1 %in% species_name) %>%
  #filter(item1 %in% gene_name) %>%
  filter(item2 %in% disease_name) %>%
  head(100)%>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = correlation, edge_width = correlation), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Top100_Relationship_by_Phi_coefficient") +
  theme_void()




# word relationship--------------------------
library(widyr)
library(igraph)
library(ggraph)

cas9_perfile.tdf <- cas9.tdf %>%
  group_by(Concept) %>%
  mutate(count_total = sum(Count)) %>%
  #summarise(count_file = n())
  #count(Concept) %>%
  arrange(desc(count_total)) %>%
  count(Concept,count_total) %>%
  arrange(desc(count_total)) 
  
table(cas9.tdf$Category)
#Choose different category for filter

abs_word_pairs <- cas9.tdf %>%
  filter(Category %in% c('Material_gene','Disease_pathologic')) %>%
  pairwise_count(Concept,PMID,sort = TRUE)

set.seed(42)
abs_word_pairs %>%
  head(100) %>%
  #filter(n < 10) %>%
  #filter(n >= 10) %>%
  #filter(n < 25) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = n, edge_width = n), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Co_occurance Concepts in abstract") +
  theme_void()

#PMI
abs_word_pairs_PMI <- cas9.tdf %>%
  filter(Category %in% c('Material_gene','Disease_pathologic')) %>%
  pairwise_pmi(Concept,PMID,sort = TRUE)


abs_word_pairs_PMI %>%
  head(100) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = pmi, edge_width = pmi), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "PMI Concepts in abstract") +
  theme_void()

abs_word_pairs_cor <- cas9.tdf %>%
  #filter(Concept %in% index) %>%
  #filter(Count >5) %>%
  pairwise_cor(Concept,PMID,Count)%>%
  arrange(desc(correlation))

abs_word_pairs_cor %>%
  #filter(correlation > 0.1) %>%
  #filter(item1 %in% c('Pancreatic Cancer')) %>%
  head(100)%>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = correlation, edge_width = correlation), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Top100_Relationship_by_Phi_coefficient") +
  theme_void()


abs_word_pairs_cosine <- cas9.tdf %>%
  #filter(Count >5) %>%
  pairwise_similarity(Concept,PMID,Count) %>%
  arrange(desc(similarity))

abs_word_pairs_cosine %>%
  head(100)%>%
  #filter(similarity > 0.5) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = similarity, edge_width = similarity), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Top100_Relationship_by_cosine_similarty") +
  theme_void()



#TFIDF --------------------------------------


tf_idf <- cas9.tdf %>%
  bind_tf_idf(PMID, Concept, Count) %>%
  arrange(desc(tf_idf))

abs_word_pairs_cosine_tfidf <- tf_idf %>%
  #filter(Count >5) %>%
  pairwise_similarity(Concept,PMID,tf_idf) %>%
  arrange(desc(similarity))

#get smtype and concept
smtype_concept <- tf_idf %>%
  count(Concept,SMtype) %>%
  group_by(Concept)


abs_word_pairs_cosine_tfidf %>%
  head(100) %>%
  #filter(similarity > 0.5) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = similarity, edge_width = similarity), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Top100_Relationship_by_cosine_similarty_TFIDF") +
  theme_void()

tf_idf_gene <- tf_idf %>%
  filter(Category == 'Material_gene')


#See the top10 tfidf in each pmid
tf_idf_top10 <- tf_idf %>%
  group_by(PMID) %>%
  filter(PMID==25775608)
  arrange(desc(PMID))%>%
  top_n(10)


#Select some PMID to show
PMID_list = sample(tf_idf$PMID,4)

tf_idf_pmid <- tf_idf %>%
  filter(PMID %in% PMID_list)
#not good view now
tf_idf_pmid %>%
  group_by(PMID) %>%
  top_n(10, tf_idf) %>%
  ungroup() %>%
  mutate(Concept = reorder(Concept, tf_idf)) %>%
  ggplot(aes(Concept, tf_idf, fill = PMID)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ PMID, scales = "free") +
  ylab("tf-idf in abstracts") +
  coord_flip()

#Shrinking and Growing words------------------------------------------
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
tmdf <- cas9.tdf %>% 
  mutate(date = as.Date(date(Date), "%Y%b%d")) %>%
  mutate(time = as.POSIXct(date, origin = "1970-01-01"))%>%
  mutate(month = round_date(date, "month"))

cas9.tdf$month = tmdf$month
#Concept count_total>100
cas9_total_concept.tdf <- cas9.tdf %>%
  count(Concept)
index <- cas9_total_concept.tdf[cas9_total_concept.tdf$n>10,]$Concept

#NEED change for PMID unique
#Change for high frequency PMID only
tmdf <- tmdf %>% 
  filter(Concept %in% index) %>%
  count(PMID,month)


papers_per_month <- tmdf %>%
  group_by(month) %>%
  summarize(month_total = n())



Concept_month_counts <- cas9.tdf %>%
  count(Concept, month) %>%
  #complete(Concept, month, Count)%>%
  complete(Concept, month, fill = list(n = 0)) %>%
  inner_join(papers_per_month, by = "month") %>%
  mutate(percent = n / month_total) %>%
  mutate(year = year(month) + yday(month) / 365) %>%
  filter(percent < 0.8)%>%
  filter(Concept %in% index) 

mod <- ~ glm(cbind(n, month_total - n) ~ year, ., family = "binomial")

slopes <- Concept_month_counts %>%
  nest(-Concept) %>%
  mutate(model = map(data, mod)) %>%
  unnest(map(model, tidy)) %>%
  filter(term == "year") %>%
  arrange(desc(estimate))

#Growing words
slopes %>%
  head(9) %>%
  inner_join(Concept_month_counts, by = "Concept") %>%
  mutate(Concept = reorder(Concept, -estimate)) %>%
  ggplot(aes(month, n / month_total, color = Concept)) +
  geom_line(show.legend = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~ Concept, scales = "free_y") +
  expand_limits(y = 0) +
  labs(x = "Year",
       y = "Percentage of titles containing this term",
       title = "9 fastest growing words",
       subtitle = "Judged by growth rate"
  )


#Shrinking words
slopes %>%
  tail(9) %>%
  inner_join(Concept_month_counts, by = "Concept") %>%
  mutate(Concept = reorder(Concept, -estimate)) %>%
  ggplot(aes(month, n / month_total, color = Concept)) +
  geom_line(show.legend = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~ Concept, scales = "free_y") +
  expand_limits(y = 0) +
  labs(x = "Year",
       y = "Percentage of titles containing this term",
       title = "9 fastest shrinking words",
       subtitle = "Judged by growth rate"
  )




# topic model--------------------------------
desc_dtm <- cas9.tdf %>%
  count(PMID, Concept, sort = TRUE) %>%
  ungroup() %>%
  cast_dtm(PMID, Concept, n)

desc_lda <- LDA(desc_dtm, k = 3, control = list(seed = 42))
tidy_lda <- tidy(desc_lda)

top_terms <- tidy_lda %>%
  group_by(topic) %>%
  top_n(10, beta) %>%
  ungroup() %>%
  arrange(topic, -beta)

top_terms %>%
  mutate(term = reorder(term, beta)) %>%
  group_by(topic, term) %>%    
  arrange(desc(beta)) %>%  
  ungroup() %>%
  mutate(term = factor(paste(term, topic, sep = "__"), 
                       levels = rev(paste(term, topic, sep = "__")))) %>%
  ggplot(aes(term, beta, fill = as.factor(topic))) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_x_discrete(labels = function(x) gsub("__.+$", "", x)) +
  labs(title = "Top 10 terms in each LDA topic",
       x = NULL, y = expression(beta)) +
  facet_wrap(~ topic, ncol = 5, scales = "free") 




