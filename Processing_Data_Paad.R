#Library Section ---------------------------

library(tidyverse)
library(ggplot2)
library(tidytext)
library(stringr)
library(scales)
library(lubridate)
library(topicmodels)

#Load Neop
paad.df <- read.delim('./paad_metamap/conceptsInFile_Neop.txt',sep = '\t')

colnames(paad.df) <- c('PMID','Date','Concept','Count')
cancer_name <- unique(paad.df$Concept)
#Load Gene
paad_gene.df <- read.delim('./paad_metamap/conceptsInFile_Gene.txt',sep = '\t')

colnames(paad_gene.df) <- c('PMID','Date','Concept','Count')

paad.tdf <- tibble::as.tibble(bind_rows(paad.df ,
                                          #filter(Concept %in% c('Pancreatic Cancer'))
                                        paad_gene.df))

#Output path
output_path <- './output/paad_metamap/'

#Shrinking and Growing words------------------------------------------
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
tmdf <- paad.tdf %>% 
  mutate(date = as.Date(date(Date), "%Y%b%d")) %>%
  mutate(time = as.POSIXct(date, origin = "1970-01-01"))%>%
  mutate(month = round_date(date, "month"))

paad.tdf$month = tmdf$month
#ZIPF_law filter count 1 concept
paad_total_concept.tdf <- paad.tdf %>%
  count(Concept)
index <- paad_total_concept.tdf[paad_total_concept.tdf$n>1,]$Concept

# word relationship--------------------------
library(widyr)
library(igraph)
library(ggraph)

abs_word_pairs <- paad.tdf %>%
  pairwise_count(Concept,PMID,sort = TRUE) %>%
  filter(item1 %in% cancer_name)
  #filter(item1 %in% c('Pancreatic Cancer'))
exclude_gene <- abs_word_pairs %>%
  filter(n==1)

exclude_gene=exclude_gene$item2

png(paste0(output_path,'ConceptRelationship.png'))
set.seed(42)
abs_word_pairs %>%
  filter(n >= 10) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = n, edge_width = n), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Disease_Gene_Relationship_by_Count") +
  theme_void()
dev.off()

#Using count for calculate corelation 
abs_word_pairs_cor <- paad.tdf %>%
  #filter(Concept %in% index) %>%
  #filter(Count >5) %>%
  pairwise_cor(Concept,PMID,Count)%>%
  arrange(desc(correlation)) %>%
  #filter(item1 %in% c('Pancreatic Cancer'))
#%>%
  filter(item1 %in% cancer_name)
#pairwise_count(Concept,PMID,sort = TRUE)
png(paste0(output_path,'ConceptRelationship_cor.png'))
abs_word_pairs_cor %>%
  filter(correlation > 0.1) %>%
  #filter(item1 %in% c('Pancreatic Cancer')) %>%
  #head(100)%>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = correlation, edge_width = correlation), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Disease_Gene_Relationship_by_Phi_coefficient") +
  theme_void()
dev.off()
# #Using count for calculate distance
# abs_word_pairs_dis <- paad.tdf %>%
#   #filter(Count >5) %>%
#   pairwise_dist(Concept,PMID,Count)%>%
#   arrange(desc(distance))
# #pairwise_count(Concept,PMID,sort = TRUE)
# 
# abs_word_pairs_dis %>%
#   filter(distance <6) %>%
#   graph_from_data_frame() %>%
#   ggraph(layout = "fr") +
#   geom_edge_link(aes(edge_alpha = distance, edge_width = distance), edge_colour = "cyan4") +
#   geom_node_point(size = 1) +
#   geom_node_text(aes(label = name), repel = TRUE, 
#                  point.padding = unit(0.2, "lines")) +
#   labs(title = "Bigrams in abstract") +
#   theme_void()


#Using count for calculate cosine similarity
abs_word_pairs_cosine <- paad.tdf %>%
  #filter(Count >5) %>%
  pairwise_similarity(Concept,PMID,Count) %>%
  arrange(desc(similarity)) %>%
  #filter(item1 %in% c('Pancreatic Cancer'))
  filter(item1 %in% cancer_name)
#pairwise_count(Concept,PMID,sort = TRUE)
png(paste0(output_path,'ConceptRelationship_cos.png'))
abs_word_pairs_cosine %>%
  head(100)%>%
  #filter(similarity > 0.5) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = similarity, edge_width = similarity), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Disease_Gene_Relationship_by_cosine_similarty") +
  theme_void()
dev.off()
#Using tfidf to calculate the relationship
tf_idf <- paad.tdf %>%
  bind_tf_idf( Concept, PMID, Count) %>%
  arrange(desc(tf))

abs_word_pairs_cosine_tfidf <- tf_idf %>%
  #filter(Count >5) %>%
  pairwise_similarity(Concept,PMID,tf_idf) %>%
  arrange(desc(similarity)) %>%
  filter(item1 %in% cancer_name)
png(paste0(output_path,'ConceptRelationship_cosine_tfidf.png'))
abs_word_pairs_cosine_tfidf %>%
  head(100) %>%
  #filter(similarity > 0.5) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = similarity, edge_width = similarity), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Disease_Gene_Relationship_by_cosine_similarty_TFIDF") +
  theme_void()
dev.off()

#NEED change for PMID unique
#Change for high frequency PMID only
tmdf <- tmdf %>% 
  filter(Concept %in% index) %>%
  count(PMID,month)


papers_per_month <- tmdf %>%
  group_by(month) %>%
  summarize(month_total = n())



Concept_month_counts <- paad.tdf %>%
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

png(paste0(output_path,'FastestGrowingWords.png'))
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
dev.off()

png(paste0(output_path,'FastestShrinkingWords.png'))
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
dev.off()

# topic model--------------------------------
desc_dtm <- paad.tdf %>%
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

png(paste0(output_path,'TopicModel_in_3.png'))
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
dev.off()



#TFIDF --------------------------------------

#Using tfidf to calculate the relationship
tf_idf <- paad.tdf %>%
  bind_tf_idf(PMID, Concept, Count) %>%
  arrange(desc(tf_idf))

abs_word_pairs_cosine_tfidf <- tf_idf %>%
  filter(Count >5) %>%
  pairwise_similarity(Concept,PMID,tf_idf)

abs_word_pairs_cosine_tfidf %>%
  filter(similarity > 0.5) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = similarity, edge_width = similarity), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Bigrams in abstract") +
  theme_void()



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
