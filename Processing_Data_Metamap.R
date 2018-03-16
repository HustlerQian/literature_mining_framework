#Library Section ---------------------------

library(tidyverse)
library(ggplot2)
library(tidytext)
library(stringr)
library(scales)
library(lubridate)
library(topicmodels)


cas9.df <- read.delim('./cas9_metamap/conceptsInFile2.txt',sep = '\t')

colnames(cas9.df) <- c('PMID','Date','Concept','Count')

cas9.tdf <- tibble::as.tibble(cas9.df)


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

# word relationship--------------------------
library(widyr)
library(igraph)
library(ggraph)

abs_word_pairs <- cas9.tdf %>%
  pairwise_count(Concept,PMID,sort = TRUE)

set.seed(42)
abs_word_pairs %>%
  #filter(n < 10) %>%
  filter(n >= 5) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = n, edge_width = n), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Bigrams in abstract") +
  theme_void()



#TFIDF --------------------------------------


tf_idf <- cas9.tdf %>%
  bind_tf_idf(PMID, Concept, Count) %>%
  arrange(desc(tf_idf))

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




