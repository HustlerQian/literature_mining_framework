#Library
library(tidyverse)
library(ggplot2)
library(tidytext)
library(stringr)
library(scales)
library(lubridate)
library(topicmodels)
#Loading data
load('./cas9.Rdata')
#Change colnames
colnames(tmdf) <- c("title"  ,  "abstract" ,"year"  ,   "month" ,   "day"    ,  "journal")  
#Output path
output_path <- './output/cas9_nonmetamap/'
#Save raw data as tmldf
tmldf <- tmdf
tmdf$line = 1:nrow(tmdf)

#Make month and date
#Change days
paperdf <- tmdf %>%
  dplyr::mutate(day = dplyr::case_when(is.na(day) ~"01",day=='NULL' ~"01", !is.na(day) ~ day))

#Change month
#Change NULL to na
index <- paperdf$month=='NULL' 
paperdf [index,'month'] = NA
paperdf = paperdf[paperdf$year!='NULL',]
#Set locales when setting the datetime
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")

paperdf <- paperdf %>%
  dplyr::mutate(month = dplyr::case_when(!(month %in% month.abb) & is.na(month) ~ "Jan", 
                                         !(month %in% month.abb) & !is.na(month) ~ month.abb[as.numeric(month)],
                                         month %in% month.abb ~ month)) %>%
  tidyr::unite(date, year, month, day, sep = "") %>%
  dplyr::mutate(date = as.Date(date, "%Y%b%d")) 

#Make time and month
tmdf <- paperdf %>% 
  mutate(time = as.POSIXct(date, origin = "1970-01-01"))%>%
  mutate(month = round_date(date, "month"))
# Numbers by journal
pdf(paste0(output_path,'Journal_paper_disturbtion.pdf'))
tmldf %>%
  group_by(journal) %>%
  summarize(papers = n_distinct(title)) %>%
  top_n(20,papers) %>%
  mutate(journal = reorder(journal, papers)) %>%     
  ggplot(aes(journal, papers)) +
  geom_col() +
  coord_flip()
dev.off()

#Title Section -------------------------------
# Get the words in title (journal with papers > 100)

pdf(paste0(output_path,'Top20UsedWordsInTitle.pdf'))
wordft <- tmdf %>%
  filter(nchar(title) > 0) %>%
  group_by(journal) %>%
  mutate(journal_total = n()) %>%
  ungroup() %>%
  filter(journal_total > 100) %>%
  unnest_tokens(word, title,drop = F) %>%
  anti_join(stop_words) %>%
  filter(str_detect(word, "[^\\d]")) %>%
  group_by(word) %>%
  mutate(word_total = n()) %>%
  ungroup() %>%
  mutate(source = 'title')

wordft %>%
  count(word, sort = TRUE) %>%
  top_n(20,n) %>%
  mutate(word = reorder(word, n)) %>%
  ggplot(aes(word, n)) +
  geom_col(show.legend = FALSE) +
  ylab("Top 20 commonly used words in titles") +
  coord_flip()
dev.off()
#Output the dataset for see
write.csv(wordft,file = paste0(output_path,'WordsFrequencyInTitle.csv'))


# Temporal Trends---------------------------
papers_per_month <- tmdf %>%
  group_by(month) %>%
  summarize(month_total = n())
# Growing words in titles
word_month_counts <- wordft %>%
  filter(word_total >= 100) %>%
  count(word, month) %>%
  complete(word, month, fill = list(n = 0)) %>%
  inner_join(papers_per_month, by = "month") %>%
  mutate(percent = n / month_total) %>%
  mutate(year = year(month) + yday(month) / 365) %>%
  filter(percent < 0.8)

mod <- ~ glm(cbind(n, month_total - n) ~ year, ., family = "binomial")

slopes <- word_month_counts %>%
  nest(-word) %>%
  mutate(model = map(data, mod)) %>%
  unnest(map(model, tidy)) %>%
  filter(term == "year") %>%
  arrange(desc(estimate))

pdf(paste0(output_path,'FastestGrowingWordsInTitle.pdf'))
slopes %>%
  head(9) %>%
  inner_join(word_month_counts, by = "word") %>%
  mutate(word = reorder(word, -estimate)) %>%
  ggplot(aes(month, n / month_total, color = word)) +
  geom_line(show.legend = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~ word, scales = "free_y") +
  expand_limits(y = 0) +
  labs(x = "Year",
       y = "Percentage of titles containing this term",
       title = "9 fastest growing words",
       subtitle = "Judged by growth rate"
  )
dev.off()

pdf(paste0(output_path,'FastestShrinkingWordsInTitle.pdf'))

# shrinking words in titles
slopes %>%
  tail(9) %>%
  inner_join(word_month_counts, by = "word") %>%
  mutate(word = reorder(word, estimate)) %>%
  ggplot(aes(month, n / month_total, color = word)) +
  geom_line(show.legend = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~ word, scales = "free_y") +
  expand_limits(y = 0) +
  labs(x = "Year",
       y = "Percentage of titles containing this term",
       title = "9 fastest shrinking words",
       subtitle = "Judged by growth rate"
  )

dev.off()


#Abstract sections----------------------------------------
pdf(paste0(output_path,'Top20UsedWordsInAbstract.pdf'))
# get words in abstracts
wordfabs <- tmdf %>%
  filter(nchar(abstract) > 0) %>%
  group_by(journal) %>%
  mutate(journal_total = n()) %>%
  ungroup() %>%
  filter(journal_total > 100) %>%
  unnest_tokens(word, abstract,drop = F) %>%
  anti_join(stop_words) %>%
  filter(str_detect(word, "[^\\d]")) %>%
  group_by(word) %>%
  mutate(word_total = n()) %>%
  ungroup() %>%
  mutate(source = 'abstract')

# plot top 20 words in abstracts
wordfabs %>%
  count(word, sort = TRUE) %>%
  top_n(20,n) %>%
  mutate(word = reorder(word, n)) %>%
  ggplot(aes(word, n)) +
  geom_col(show.legend = FALSE) +
  ylab("Top 20 commonly used words in abstracts") +
  coord_flip()
#Output the dataset for see
write.csv(wordfabs,file = paste0(output_path,'WordsFrequencyInAbstract.csv'))
dev.off()



# Growing words in abstracts
word_month_counts <- wordfabs %>%
  filter(word_total >= 100) %>%
  count(word, month) %>%
  complete(word, month, fill = list(n = 0)) %>%
  inner_join(papers_per_month, by = "month") %>%
  mutate(percent = n / month_total) %>%
  mutate(year = year(month) + yday(month) / 365) %>%
  filter(percent < 0.8)

mod <- ~ glm(cbind(n, month_total - n) ~ year, ., family = "binomial")

pdf(paste0(output_path,'FastestGrowingWordsInAbstract.pdf'))
slopes <- word_month_counts %>%
  nest(-word) %>%
  mutate(model = map(data, mod)) %>%
  unnest(map(model, tidy)) %>%
  filter(term == "year") %>%
  arrange(desc(estimate))

slopes %>%
  head(9) %>%
  inner_join(word_month_counts, by = "word") %>%
  mutate(word = reorder(word, -estimate)) %>%
  ggplot(aes(month, n / month_total, color = word)) +
  geom_line(show.legend = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~ word, scales = "free_y") +
  expand_limits(y = 0) +
  labs(x = "Year",
       y = "Percentage of abstracts containing this term",
       title = "9 fastest growing words ",
       subtitle = "Judged by growth rate")
dev.off()
# shrinking words in abstracts
pdf(paste0(output_path,'FastestShrinkingWordsInAbstract.pdf'))
slopes %>%
  tail(9) %>%
  inner_join(word_month_counts, by = "word") %>%
  mutate(word = reorder(word, estimate)) %>%
  ggplot(aes(month, n / month_total, color = word)) +
  geom_line(show.legend = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~ word, scales = "free_y") +
  expand_limits(y = 0) +
  labs(x = "Year",
       y = "Percentage of abstracts containing this term",
       title = "9 fastest shrinking words",
       subtitle = "Judged by growth rate")
dev.off()

#Word relationship -----------------------------------------------
library(widyr)
library(igraph)
library(ggraph)

#Title word relationship
title_word_pairs <- wordft %>%
  pairwise_count(word,line,sort = TRUE)

#Filter the correct times
set.seed(42)
pdf(paste0(output_path,'WordsRelationshipInTitle.pdf'))
title_word_pairs %>%
  filter(n >= 20) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = n, edge_width = n), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Bigrams in title") +
  theme_void()
dev.off()

#Abstract word relationship
abs_word_pairs <- wordfabs %>%
pairwise_count(word,line,sort = TRUE)
pdf(paste0(output_path,'WordsRelationshipInAbstract.pdf'))
set.seed(42)
abs_word_pairs %>%
  filter(n >= 100) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = n, edge_width = n), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Bigrams in abstract") +
  theme_void()
dev.off()

# topic model--------------------------------
desc_dtm <- wordfabs %>%
  count(line, word, sort = TRUE) %>%
  ungroup() %>%
  cast_dtm(line, word, n)

desc_lda <- LDA(desc_dtm, k = 5, control = list(seed = 42))
tidy_lda <- tidy(desc_lda)

top_terms <- tidy_lda %>%
  group_by(topic) %>%
  top_n(10, beta) %>%
  ungroup() %>%
  arrange(topic, -beta)

#Topic models
pdf(paste0(output_path,'TopicModel_in_5.pdf'))
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

# Journal accent---------------------------

words_by_journal <- wordfabs %>%
  count(journal, word, sort = TRUE) %>%
  ungroup()

tf_idf <- words_by_journal %>%
  bind_tf_idf(word, journal, n) %>%
  arrange(desc(tf_idf))

tf_idf %>%
  group_by(journal) %>%
  top_n(10, tf_idf) %>%
  ungroup() %>%
  mutate(word = reorder(word, tf_idf)) %>%
  ggplot(aes(word, tf_idf, fill = journal)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ journal, scales = "free") +
  ylab("tf-idf in abstracts") +
  coord_flip()

#Sentimental Analysis----------------------

contributions <- wordft %>%
  inner_join(get_sentiments("afinn"), by = "word") %>%
  group_by(word) %>%
  summarize(occurences = n(),
            contribution = sum(score))

contributions %>%
  top_n(25, abs(contribution)) %>%
  mutate(word = reorder(word, contribution)) %>%
  ggplot(aes(word, contribution, fill = contribution > 0)) +
  geom_col(show.legend = FALSE) +
  coord_flip()

#Abstract
contributions <- wordfabs %>%
  inner_join(get_sentiments("afinn"), by = "word") %>%
  group_by(word) %>%
  summarize(occurences = n(),
            contribution = sum(score))

contributions %>%
  top_n(25, abs(contribution)) %>%
  mutate(word = reorder(word, contribution)) %>%
  ggplot(aes(word, contribution, fill = contribution > 0)) +
  geom_col(show.legend = FALSE) +
  coord_flip()


