#Library Section ---------------------------

library(tidyverse)
library(ggplot2)
library(tidytext)
library(stringr)
library(scales)
library(lubridate)
library(topicmodels)
library(widyr)
library(igraph)
library(ggraph)

options(stringsAsFactors = F)

#Load PAAD-----------------------------------
paad.df <- read.delim('./paad_metamap/conceptsInFile_Neop.txt',sep = '\t')

colnames(paad.df) <- c('PMID','Date','Concept','Count')
cancer_name <- unique(paad.df$Concept)
#Load Gene
paad_gene.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/paad/conceptsInFile_Gene.txt',sep = '\t')

colnames(paad_gene.df) <- c('PMID','Date','Concept','Count')

#Load gene2 
paad_gene2.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/paad/conceptsInFile_Gene2.txt',sep = '\t')

colnames(paad_gene2.df) <- c('PMID','Date','Concept','Count')

paad_gene.df <- tibble::as.tibble(bind_rows(paad_gene.df,
                                            paad_gene2.df))
paad.tdf <- tibble::as.tibble(bind_rows(paad.df ,
                                        #filter(Concept %in% c('Pancreatic Cancer'))
                                        paad_gene.df))

#Load chol-----------------------------------
chol.df <- read.delim('./chol_metamap/conceptsInFile_Neop.txt',sep = '\t')

colnames(chol.df) <- c('PMID','Date','Concept','Count')
cancer_name <- unique(chol.df$Concept)
#Load Gene
chol_gene.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/chol/conceptsInFile_Gene.txt',sep = '\t')

colnames(chol_gene.df) <- c('PMID','Date','Concept','Count')
#Load gene2 
chol_gene2.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/chol/conceptsInFile_Gene2.txt',sep = '\t')

colnames(chol_gene2.df) <- c('PMID','Date','Concept','Count')

chol_gene.df <- tibble::as.tibble(bind_rows(chol_gene.df,
                                            chol_gene2.df))

chol.tdf <- tibble::as.tibble(bind_rows(chol.df ,
                                        #filter(Concept %in% c('Pancreatic Cancer'))
                                        chol_gene.df))
#Load esca-----------------------------------
esca.df <- read.delim('./esca_metamap/conceptsInFile_Neop.txt',sep = '\t')

colnames(esca.df) <- c('PMID','Date','Concept','Count')
cancer_name <- unique(esca.df$Concept)
#Load Gene
esca_gene.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/esca/conceptsInFile_Gene.txt',sep = '\t')

colnames(esca_gene.df) <- c('PMID','Date','Concept','Count')
#Load gene2 
esca_gene2.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/esca/conceptsInFile_Gene2.txt',sep = '\t')

colnames(esca_gene2.df) <- c('PMID','Date','Concept','Count')

esca_gene.df <- tibble::as.tibble(bind_rows(esca_gene.df,
                                            esca_gene2.df))
esca.tdf <- tibble::as.tibble(bind_rows(esca.df ,
                                        #filter(Concept %in% c('Pancreatic Cancer'))
                                        esca_gene.df))
#Load coca-----------------------------------
coca.df <- read.delim('./coca_metamap/conceptsInFile_Neop.txt',sep = '\t')

colnames(coca.df) <- c('PMID','Date','Concept','Count')
cancer_name <- unique(coca.df$Concept)
#Load Gene
coca_gene.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/coca/conceptsInFile_Gene.txt',sep = '\t')

colnames(coca_gene.df) <- c('PMID','Date','Concept','Count')
#Load gene2 
coca_gene2.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/coca/conceptsInFile_Gene2.txt',sep = '\t')

colnames(coca_gene2.df) <- c('PMID','Date','Concept','Count')

coca_gene.df <- tibble::as.tibble(bind_rows(coca_gene.df,
                                            coca_gene2.df))

coca.tdf <- tibble::as.tibble(bind_rows(coca.df ,
                                        #filter(Concept %in% c('Pancreatic Cancer'))
                                        coca_gene.df))

#Load stad-----------------------------------
stad.df <- read.delim('./stad_metamap/conceptsInFile_Neop.txt',sep = '\t')

colnames(stad.df) <- c('PMID','Date','Concept','Count')
cancer_name <- unique(stad.df$Concept)
#Load Gene
stad_gene.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/stad/conceptsInFile_Gene.txt',sep = '\t')

colnames(stad_gene.df) <- c('PMID','Date','Concept','Count')
#Load gene2 
stad_gene2.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/stad/conceptsInFile_Gene2.txt',sep = '\t')

colnames(stad_gene2.df) <- c('PMID','Date','Concept','Count')

stad_gene.df <- tibble::as.tibble(bind_rows(stad_gene.df,
                                            stad_gene2.df))

stad.tdf <- tibble::as.tibble(bind_rows(stad.df ,
                                        #filter(Concept %in% c('Pancreatic Cancer'))
                                        stad_gene.df))
#Load lihc-----------------------------------
lihc.df <- read.delim('./lihc_metamap/conceptsInFile_Neop.txt',sep = '\t')

colnames(lihc.df) <- c('PMID','Date','Concept','Count')
cancer_name <- unique(lihc.df$Concept)
#Load Gene
lihc_gene.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/lihc/conceptsInFile_Gene.txt',sep = '\t')

colnames(lihc_gene.df) <- c('PMID','Date','Concept','Count')
#Load gene2 
lihc_gene2.df <- read.delim('D:/wyq/Gi_Gene_Timeline/res/lihc/conceptsInFile_Gene2.txt',sep = '\t')

colnames(lihc_gene2.df) <- c('PMID','Date','Concept','Count')

lihc_gene.df <- tibble::as.tibble(bind_rows(lihc_gene.df,
                                            lihc_gene2.df))

lihc.tdf <- tibble::as.tibble(bind_rows(lihc.df ,
                                        #filter(Concept %in% c('Pancreatic Cancer'))
                                        lihc_gene.df))


#Combine together 
GI.tdf <- bind_rows(stad.tdf,paad.tdf,lihc.tdf,esca.tdf,coca.tdf,chol.tdf)



output_path <- './output/GI_metamap/'

#Output　the concept freq in every file
write.csv(GI.tdf,file = paste0(output_path,'GI_all_concept_frequency.csv'))

#Exclude gene concept
gene_exclude <- c('TAT','HCCS','BTC',
                  'CCS','HPD','FAP','MUT','CAT','PIGS','ITCH','SDS','GBA','PCCA','KIT','SHE',
                  'LARGE','NIN','GC','BID','VIP','RAN','POLE','AGA','GIF','TH','EDA','NNT','VIT',
                  'LCT','HGS','CIC','DST','NRK','TTN','NTM','MBP','EFS ','INA','F12','HDC','HPN',
                  'PLN','STS','TOX','PC','INS','MAK','HGD','LGD','ABO','AGT','F2','MIP','CGA','NPS',
                  'HAL','MIA','F3','PTS','F9','F10','F7','RDX',
                  'F11','LPL','F8','KEL','NEB','SELL','MLN','PAM','GCG','TUB','HPR','PSD','BLM','GNE','PTN')
GI.tdf <- GI.tdf %>%
  filter(!Concept %in% gene_exclude)

#use corelation----------------------
#png(paste0(output_path,'ConceptRelationship.png'))
pdf(paste0(output_path,'ConceptRelationship.pdf'))
abs_word_pairs <- GI.tdf %>%
  pairwise_count(Concept,PMID,sort = TRUE) %>%
  filter(item1 %in% cancer_name) %>%
  filter(!item2 %in% cancer_name) %>%
  rename(Cancer=item1,Gene=item2)

write.csv(abs_word_pairs,file = paste0(output_path,'GI_disease_gene_paircount.csv'),row.names = F)
#Use spread to see the gene type
abs_word_pairs_sp <- abs_word_pairs %>% spread(Cancer,n) %>% 
#replace na to zero : https://stackoverflow.com/questions/45576805/how-to-replace-all-na-in-a-dataframe-using-tidyrreplace-na
          replace(is.na(.),0) 
#Show colnames
colnames(abs_word_pairs_sp)
abs_word_pairs_sp <- abs_word_pairs_sp %>%
          mutate(total = `Bile duct cancer`+`Colorectal cancer` + `Esophageal cancer`+ `Liver cancer`+`Pancreatic Cancer` + `Stomach cancer`) %>%
          arrange(desc(total))
write.csv(abs_word_pairs_sp,file = paste0(output_path,'GI_disease_gene_paircount_sp.csv'),row.names = F)

abs_word_pairs %>%
  head(100) %>%
  #filter(n >= 50) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = n, edge_width = n), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE,size=5,
                 point.padding = unit(0.2, "lines")) +
  theme(legend.text = element_text(size=20))+
  #labs(title = "Disease_Gene_Relationship_by_Count") +
  theme_void()
#dev.off()
#Get top 100
# Top100 <- abs_word_pairs %>%
#   head(100)
# 
# #Make nodes datasets
# nodes <- data.frame(name=c(unique(Top100$Cancer),unique(Top100$Gene)),
#                     group=c(rep('Cancer',length(unique(Top100$Cancer))),
#                             rep('Gene',length(unique(Top100$Gene)))
#                     )
# )
# 
# #Make graph dataframe
# net <- graph_from_data_frame(d = Top100,vertices = nodes, directed=T)
# 
# net %>%
#   ggraph(layout = "fr") +
#   geom_edge_link(aes(edge_alpha = n, edge_width = n), edge_colour = "cyan4") +
#   geom_node_point(aes(col = group),size = 2 ) +
#   geom_node_text(aes(label = name), repel = TRUE, 
#                  point.padding = unit(0.2, "lines")) +
#   labs(title = "Disease_Gene_Relationship_by_Count") +
#   theme_void()

#Use PMI -----------------------------
abs_word_pairs_PMI <- GI.tdf %>%
  pairwise_pmi(Concept,PMID,sort = TRUE) %>%
  filter(item1 %in% cancer_name) %>%
  filter(!item2 %in% cancer_name) %>%
  rename(Cancer=item1,Gene=item2)
write.csv(abs_word_pairs_PMI,file = paste0(output_path,'GI_disease_gene_pairPMI.csv'),row.names = F)

abs_word_pairs_PMI_sp <- abs_word_pairs_PMI %>% spread(Cancer,PMI) %>% 
  #replace na to zero : https://stackoverflow.com/questions/45576805/how-to-replace-all-na-in-a-dataframe-using-tidyrreplace-na
  replace(is.na(.),0) 
#Show colnames
colnames(abs_word_pairs_PMI_sp)
abs_word_pairs_PMI_sp <- abs_word_pairs_PMI_sp %>%
  mutate(total = `Bile duct cancer`+`Colorectal cancer` + `Esophageal cancer`+ `Liver cancer`+`Pancreatic Cancer` + `Stomach cancer`) %>%
  arrange(desc(total))
write.csv(abs_word_pairs_PMI_sp,file = paste0(output_path,'GI_disease_gene_PMI_sp.csv'),row.names = F)

abs_word_pairs_PMI %>%
  #filter(n >= 50) %>%
  head(100) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = pmi, edge_width = pmi), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, size=6,
                 point.padding = unit(0.2, "lines")) +
  #labs(title = "Disease_Gene_Relationship_by_PMI") +
  theme_void()


#Use corelation--------------------
#Using count for calculate corelation 
abs_word_pairs_cor <- GI.tdf %>%
  pairwise_cor(Concept,PMID,Count)%>%
  arrange(desc(correlation)) %>%
  filter(item1 %in% cancer_name) %>%
  filter(!item2 %in% cancer_name) %>%
  rename(Cancer=item1,Gene=item2)
write.csv(abs_word_pairs_cor,file = paste0(output_path,'GI_disease_gene_paircorelation.csv'),row.names = F)

abs_word_pairs_cor_sp <- abs_word_pairs_cor %>% spread(Cancer,correlation) %>% 
  #replace na to zero : https://stackoverflow.com/questions/45576805/how-to-replace-all-na-in-a-dataframe-using-tidyrreplace-na
  replace(is.na(.),0) 
#Show colnames
colnames(abs_word_pairs_cor_sp)
abs_word_pairs_cor_sp <- abs_word_pairs_cor_sp %>%
  mutate(total = `Bile duct cancer`+`Colorectal cancer` + `Esophageal cancer`+ `Liver cancer`+`Pancreatic Cancer` + `Stomach cancer`) %>%
  arrange(desc(total))
write.csv(abs_word_pairs_cor_sp,file = paste0(output_path,'GI_disease_gene_cor_sp.csv'),row.names = F)


#png(paste0(output_path,'ConceptRelationship_cor.png'))
abs_word_pairs_cor %>%
  #filter(!item2 %in% cancer_name)%>%
  #filter(correlation > 0.1) %>%
  #filter(item1 %in% c('Pancreatic Cancer')) %>%
  head(100)%>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = correlation, edge_width = correlation), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, size=5,
                 point.padding = unit(0.2, "lines")) +
  #labs(title = "Disease_Gene_Relationship_by_Phi_coefficient") +
  theme_void()

#dev.off()

#Use cosine---------------------
#png(paste0(output_path,'ConceptRelationship_cos.png'))
abs_word_pairs_cosine <- GI.tdf %>%
  pairwise_similarity(Concept,PMID,Count) %>%
  arrange(desc(similarity)) %>%
  filter(item1 %in% cancer_name)%>%
  filter(!item2 %in% cancer_name)%>%
  rename(Cancer=item1,Gene=item2)
write.csv(abs_word_pairs_cosine,file = paste0(output_path,'GI_disease_gene_paircosine.csv'),row.names = F)
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
#dev.off()

tf_idf <- GI.tdf %>%
  bind_tf_idf( Concept, PMID, Count) %>%
  arrange(desc(tf_idf))
#Save tfidf result for show
write.csv(tf_idf,file='GI_tfidf.csv',row.names = F)

abs_word_pairs_cosine_tfidf <- tf_idf %>%
  #filter(Count >5) %>%
  pairwise_similarity(Concept,PMID,tf_idf) %>%
  arrange(desc(similarity)) %>%
  filter(item1 %in% cancer_name)%>%
  filter(!item2 %in% cancer_name)%>%
  rename(Cancer=item1,Gene=item2)
write.csv(abs_word_pairs_cosine_tfidf,file = paste0(output_path,'GI_disease_gene_paircosine_tfidf.csv'),row.names = F)

#png(paste0(output_path,'ConceptRelationship_cosine_tfidf.png'))
abs_word_pairs_cosine_tfidf %>%
  head(100) %>%
  #filter(similarity > 0.5) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = similarity, edge_width = similarity), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, size=5,
                 point.padding = unit(0.2, "lines")) +
  #labs(title = "Disease_Gene_Relationship_by_cosine_similarty_TFIDF") +
  theme_void()
dev.off()


#End-------------------------------------------
png(paste0(output_path,'ConceptRelationship_cosine_tf.png'))
abs_word_pairs_cosine_tf <- tf_idf %>%
  #filter(Count >5) %>%
  pairwise_similarity(Concept,PMID,tf) %>%
  arrange(desc(similarity)) %>%
  filter(item1 %in% cancer_name)

abs_word_pairs_cosine_tf %>%
  head(100) %>%
  #filter(similarity > 0.5) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = similarity, edge_width = similarity), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Disease_Gene_Relationship_by_cosine_similarty_TF") +
  theme_void()
dev.off()
png(paste0(output_path,'ConceptRelationship_cosine_idf.png'))
abs_word_pairs_cosine_idf <- tf_idf %>%
  #filter(Count >5) %>%
  pairwise_similarity(Concept,PMID,idf) %>%
  arrange(desc(similarity)) %>%
  filter(item1 %in% cancer_name)

abs_word_pairs_cosine_idf %>%
  head(100) %>%
  #filter(similarity > 0.5) %>%
  graph_from_data_frame() %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(edge_alpha = similarity, edge_width = similarity), edge_colour = "cyan4") +
  geom_node_point(size = 1) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 point.padding = unit(0.2, "lines")) +
  labs(title = "Disease_Gene_Relationship_by_cosine_similarty_IDF") +
  theme_void()
dev.off()

