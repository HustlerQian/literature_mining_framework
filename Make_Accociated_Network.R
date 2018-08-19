library(tidyverse)
library(tidytext)
library(widyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(scales)
library(Cairo)
setwd('D:/GitHub/literature_mining_framework/')
output_path <- "./output/GI_metamap/"
load('./GI.RData')

write.csv(GI.tdf[,-2],row.names = F,
          file = paste0(output_path,'GI_all_concept_frequency.csv'))

#Make cooccurence table----------------------------------------------
abs_word_pairs <- GI.tdf %>%
  pairwise_count(Concept,PMID,sort = TRUE) %>%
  filter(item1 %in% cancer_name) %>%
  filter(!item2 %in% cancer_name) %>%
  rename(Cancer=item1,Gene=item2)

abs_word_pairs_sp <- abs_word_pairs %>% spread(Cancer,n) %>% 
  #replace na to zero : https://stackoverflow.com/questions/45576805/how-to-replace-all-na-in-a-dataframe-using-tidyrreplace-na
  replace(is.na(.),0) 
colnames(abs_word_pairs_sp)
abs_word_pairs_sp <- abs_word_pairs_sp %>%
  mutate(total = `Bile duct cancer`+`Colorectal cancer` + `Esophageal cancer`+ `Liver cancer`+`Pancreatic Cancer` + `Stomach cancer`) %>%
  arrange(desc(total))
write.csv(abs_word_pairs_sp,file = paste0(output_path,'GI_disease_gene_paircount_sp.csv'),row.names = F)


#For cooccurence network----------------------------------------------------------




abs_word_pairs <- GI.tdf %>%
  pairwise_count(Concept,PMID,sort = TRUE) %>%
  filter(item1 %in% cancer_name) %>%
  filter(!item2 %in% cancer_name) %>%
  rename(Cancer=item1,Gene=item2,Cooccur=n) %>%
  head(100)

#Make node 
cancer=unique(abs_word_pairs$Cancer)
gene=unique(abs_word_pairs$Gene)
count_node=data.frame(name=c(cancer,gene),
                     group=c(rep('cancer',length(cancer)),rep('gene',length(gene))),
                     size=10
)
count_node$id=0:(nrow(count_node)-1)
count_node$label=count_node$name
count_node_dic=count_node
rownames(count_node_dic)=count_node_dic$name
#Make edge
count_edge <- abs_word_pairs
count_edge$from <- count_node_dic[abs_word_pairs$Cancer,'id']
count_edge$to <- count_node_dic[abs_word_pairs$Gene,'id']


mygraph <- graph_from_data_frame( count_edge, vertices=count_node)
mygraph_tidy <- mygraph %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(gene=!node_is_source())%>%
  mutate(cancer=node_is_source())%>%
  activate(edges)
mygraph_tidy %>%
  ggraph(layout = 'fr') +
  geom_edge_link(aes(edge_alpha = Cooccur, edge_width = Cooccur), edge_colour = "green") +
  geom_node_point(aes(filter=cancer),color='#F45B24',size =5)+
  geom_node_point(aes(filter=gene),color='#547B94',size =2)+
  geom_node_text(aes(filter=gene,label = name),color='#547B94',repel = TRUE,
                 point.padding = unit(0.2, "lines"))+
  geom_node_text(aes(filter=cancer,label = name),color='#F45B24',repel = TRUE,
                 point.padding = unit(0.2, "lines"))+
  theme_graph()

#PHI----------------------------------------------------
abs_word_pairs_phi <- GI.tdf %>%
  pairwise_cor(Concept,PMID,Count)%>%
  arrange(desc(correlation)) %>%
  filter(item1 %in% cancer_name) %>%
  filter(!item2 %in% cancer_name) %>%
  rename(Cancer=item1,Gene=item2,Phi=correlation)

abs_word_pairs_phi_sp <- abs_word_pairs_phi %>% spread(Cancer,Phi) %>% 
  #replace na to zero : https://stackoverflow.com/questions/45576805/how-to-replace-all-na-in-a-dataframe-using-tidyrreplace-na
  replace(is.na(.),0) 
colnames(abs_word_pairs_phi_sp)
abs_word_pairs_phi_sp <- abs_word_pairs_phi_sp %>%
  mutate(total = `Bile duct cancer`+`Colorectal cancer` + `Esophageal cancer`+ `Liver cancer`+`Pancreatic Cancer` + `Stomach cancer`) %>%
  arrange(desc(total))
write.csv(abs_word_pairs_phi_sp,file = paste0(output_path,'GI_disease_gene_pairPHI_sp.csv'),row.names = F)

abs_word_pairs_phi <- GI.tdf %>%
  pairwise_cor(Concept,PMID,Count)%>%
  arrange(desc(correlation)) %>%
  filter(item1 %in% cancer_name) %>%
  filter(!item2 %in% cancer_name) %>%
  rename(Cancer=item1,Gene=item2,Phi=correlation) %>%
  head(100)

#Make node 
cancer=unique(abs_word_pairs_phi$Cancer)
gene=unique(abs_word_pairs_phi$Gene)
phi_node=data.frame(name=c(cancer,gene),
                    group=c(rep('cancer',length(cancer)),rep('gene',length(gene))),
                    size=10
)
phi_node$id=0:(nrow(phi_node)-1)
phi_node$label=phi_node$name
phi_node_dic=phi_node
rownames(phi_node_dic)=phi_node_dic$name
#Make edge
phi_edge <- abs_word_pairs_phi
phi_edge$from <- phi_node_dic[abs_word_pairs_phi$Cancer,'id']
phi_edge$to <- phi_node_dic[abs_word_pairs_phi$Gene,'id']


mygraph <- graph_from_data_frame( phi_edge, vertices=phi_node)
mygraph_tidy <- mygraph %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(gene=!node_is_source())%>%
  mutate(cancer=node_is_source())%>%
  activate(edges)
mygraph_tidy %>%
  ggraph(layout = 'fr') +
  geom_edge_link(aes(edge_alpha = Phi, edge_width = Phi), edge_colour = "green") +
  geom_node_point(aes(filter=cancer),color='#F45B24',size =5)+
  geom_node_point(aes(filter=gene),color='#547B94',size =2)+
  geom_node_text(aes(filter=gene,label = name),color='#547B94',repel = TRUE,
                 point.padding = unit(0.2, "lines"))+
  geom_node_text(aes(filter=cancer,label = name),color='#F45B24',repel = TRUE,
                 point.padding = unit(4, "lines"))+
  theme_graph()


#PMI -----------------------------------------------
abs_word_pairs_PMI <- GI.tdf %>%
  pairwise_pmi(Concept,PMID,sort = TRUE) %>%
  filter(item1 %in% cancer_name) %>%
  filter(!item2 %in% cancer_name) %>%
  rename(Cancer=item1,Gene=item2,PMI=pmi) 

abs_word_pairs_PMI_sp <- abs_word_pairs_PMI %>% spread(Cancer,PMI) %>% 
  #replace na to zero : https://stackoverflow.com/questions/45576805/how-to-replace-all-na-in-a-dataframe-using-tidyrreplace-na
  replace(is.na(.),0) 
colnames(abs_word_pairs_PMI_sp)
abs_word_pairs_PMI_sp <- abs_word_pairs_PMI_sp %>%
  mutate(total = `Bile duct cancer`+`Colorectal cancer` + `Esophageal cancer`+ `Liver cancer`+`Pancreatic Cancer` + `Stomach cancer`) %>%
  arrange(desc(total))
write.csv(abs_word_pairs_PMI_sp,file = paste0(output_path,'GI_disease_gene_pairPMI_sp.csv'),row.names = F)


abs_word_pairs_pmi <-  abs_word_pairs_PMI %>%
    head(100)

#Make node 
cancer=unique(abs_word_pairs_pmi$Cancer)
gene=unique(abs_word_pairs_pmi$Gene)
pmi_node=data.frame(name=c(cancer,gene),
                    group=c(rep('cancer',length(cancer)),rep('gene',length(gene))),
                    size=10
)
pmi_node$id=0:(nrow(pmi_node)-1)
pmi_node$label=pmi_node$name
pmi_node_dic=pmi_node
rownames(pmi_node_dic)=pmi_node_dic$name
#Make edge
pmi_edge <- abs_word_pairs_pmi
pmi_edge$from <- pmi_node_dic[abs_word_pairs_pmi$Cancer,'id']
pmi_edge$to <- pmi_node_dic[abs_word_pairs_pmi$Gene,'id']


mygraph <- graph_from_data_frame( pmi_edge, vertices=pmi_node)
mygraph_tidy <- mygraph %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(gene=!node_is_source())%>%
  mutate(cancer=node_is_source())%>%
  activate(edges)
mygraph_tidy %>%
  ggraph(layout = 'fr') +
  geom_edge_link(aes(edge_alpha = PMI, edge_width = PMI), edge_colour = "green") +
  geom_node_point(aes(filter=cancer),color='#F45B24',size =5)+
  geom_node_point(aes(filter=gene),color='#547B94',size =2)+
  geom_node_text(aes(filter=gene,label = name),color='#547B94',repel = TRUE,
                 point.padding = unit(0.2, "lines"))+
  geom_node_text(aes(filter=cancer,label = name),color='#F45B24',repel = TRUE,
                 point.padding = unit(1, "lines"))+
  theme_graph()

PMI_100 <- head(abs_word_pairs_PMI,n=100)
#Boxplot
ggplot(PMI_100,aes(x=Cancer,y=pmi))+
  geom_boxplot()
#PieChart
PMI_100_stat <- PMI_100 %>% 
  group_by(Cancer) %>%
  summarise(value =n()) %>%
  arrange(desc(value))

#Spread for heatmap
abs_word_pairs_PMI_sp <- abs_word_pairs_PMI %>% 
  spread(Cancer,pmi) %>%
  replace(is.na(.),0)
rownames(abs_word_pairs_PMI_sp)=abs_word_pairs_PMI_sp$Gene
pheatmap::pheatmap(abs_word_pairs_PMI_sp[,-1],cluster_rows = F,cluster_cols = F)


#Cosine Similarity----------------------------------------------
abs_word_pairs_cosine <- GI.tdf %>%
  pairwise_similarity(Concept,PMID,Count) %>%
  arrange(desc(similarity)) %>%
  filter(item1 %in% cancer_name)%>%
  filter(!item2 %in% cancer_name)%>%
  rename(Cancer=item1,Gene=item2,Cosine=similarity)

#Note tfidf trans with cosine
abs_word_pairs_cosine_sp <- abs_word_pairs_cosine %>% spread(Cancer,Cosine) %>% 
  #replace na to zero : https://stackoverflow.com/questions/45576805/how-to-replace-all-na-in-a-dataframe-using-tidyrreplace-na
  replace(is.na(.),0) 
colnames(abs_word_pairs_cosine_sp)
abs_word_pairs_cosine_sp <- abs_word_pairs_cosine_sp %>%
  mutate(total = `Bile duct cancer`+`Colorectal cancer` + `Esophageal cancer`+ `Liver cancer`+`Pancreatic Cancer` + `Stomach cancer`) %>%
  arrange(desc(total))
write.csv(abs_word_pairs_cosine_sp,file = paste0(output_path,'GI_disease_gene_pairCosineTFIDF_sp.csv'),row.names = F)


abs_word_pairs_cosine <- GI.tdf %>%
  pairwise_similarity(Concept,PMID,Count) %>%
  arrange(desc(similarity)) %>%
  filter(item1 %in% cancer_name)%>%
  filter(!item2 %in% cancer_name)%>%
  rename(Cancer=item1,Gene=item2,Cosine=similarity) %>%
  head(100)

#Add tfidf-------------------------------------------------
library(tidytext)
tf_idf <- GI.tdf %>%
  bind_tf_idf( Concept, PMID, Count) %>%
  arrange(desc(tf_idf))

abs_word_pairs_cosine_tfidf <- tf_idf %>%
  pairwise_similarity(Concept,PMID,tf_idf) %>%
  arrange(desc(similarity)) %>%
  filter(item1 %in% cancer_name)%>%
  filter(!item2 %in% cancer_name)%>%
  rename(Cancer=item1,Gene=item2,Cosine=similarity) 
#%>%
  #head(100)


write.csv(abs_word_pairs_cosine_tfidf,file = paste0(output_path,'GI_disease_gene_paircosine_tfidf.csv'),row.names = F)

abs_word_pairs_cosine_tfidf_sp <- abs_word_pairs_cosine_tfidf %>% spread(Cancer,Cosine) %>% 
  #replace na to zero : https://stackoverflow.com/questions/45576805/how-to-replace-all-na-in-a-dataframe-using-tidyrreplace-na
  replace(is.na(.),0) 
colnames(abs_word_pairs_cosine_tfidf_sp)
abs_word_pairs_cosine_tfidf_sp <- abs_word_pairs_cosine_tfidf_sp %>%
  mutate(total = `Bile duct cancer`+`Colorectal cancer` + `Esophageal cancer`+ `Liver cancer`+`Pancreatic Cancer` + `Stomach cancer`) %>%
  arrange(desc(total))
write.csv(abs_word_pairs_cosine_tfidf_sp,file = paste0(output_path,'GI_disease_gene_pairCosine_sp.csv'),row.names = F)

#Make node 
cancer=unique(abs_word_pairs_cosine$Cancer)
gene=unique(abs_word_pairs_cosine$Gene)
cosine_node=data.frame(name=c(cancer,gene),
                       group=c(rep('cancer',length(cancer)),rep('gene',length(gene))),
                       size=10
)
cosine_node$id=0:(nrow(cosine_node)-1)
cosine_node$label=cosine_node$name
cosine_node_dic=cosine_node
rownames(cosine_node_dic)=cosine_node_dic$name
#Make edge
cosine_edge <- abs_word_pairs_cosine
cosine_edge$from <- cosine_node_dic[abs_word_pairs_cosine$Cancer,'id']
cosine_edge$to <- cosine_node_dic[abs_word_pairs_cosine$Gene,'id']


mygraph <- graph_from_data_frame( cosine_edge, vertices=cosine_node)
mygraph_tidy <- mygraph %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(gene=!node_is_source())%>%
  mutate(cancer=node_is_source())%>%
  activate(edges)
mygraph_tidy %>%
  ggraph(layout = 'fr') +
  geom_edge_link(aes(edge_alpha = Cosine, edge_width = Cosine), edge_colour = "green") +
  geom_node_point(aes(filter=cancer),color='#F45B24',size =5)+
  geom_node_point(aes(filter=gene),color='#547B94',size =2)+
  geom_node_text(aes(filter=gene,label = name),color='#547B94',repel = TRUE,
                 point.padding = unit(0.2, "lines"))+
  geom_node_text(aes(filter=cancer,label = name),color='#F45B24',repel = TRUE,
                 point.padding = unit(1, "lines"))+
  theme_graph()

#pdf('Three_methods_spearman_coeffient_by_cancer_types.pdf')

# pdf(file = paste0(normalizePath('D:/wyq/消化道肿瘤热点基因文献挖掘/pics/'),
#            'Spearman_Correlationship_four_methods.pdf'),
#     paper = 'a4r')
#Correlationship plot for 4 method and 6 cancers.-----------------------
#Change the color palette
#col <- colorRampPalette(c("blue", "white", "red"))

col <- colorRampPalette(c("#7F0000", "red", 
                          "white",
                         "blue","#00007F"))
#Change the par margin form of layout
par(mfrow = c(2,2))


#Spearman PHI correlationship for by different cancer type-------------------------
library(corrplot)

#View(abs_word_pairs_phi_sp)
cor.m <- abs_word_pairs_phi_sp[,-1]
colnames(cor.m) <- c("Bile duct" ,
                     "Colorectal" ,
                     "Esophageal" ,
                     "Liver" , 
                     "Pancreatic" ,
                     "Stomach")
Phi_cor_spearman <- signif(cor(cor.m,method = "spearman"),3)

corrplot(Phi_cor_spearman,
         method = "color",
         type = "upper",
         col = col(50),
         tl.pos = "d",
         mar=c(0,0,1.5,1),
         title = 'A. Phi spearman coeffient',
         cl.length = 11
         )

corrplot(Phi_cor_spearman, 
         add = TRUE, 
         type = "lower", 
         method = "number", 
         diag = FALSE, 
         col = "black",
         tl.pos = "n", 
         mar=c(0,0,1.5,1),
         cl.pos = "n")
#Old style for backup---------------------
# dev.off()
# text(1,1,'Phi spearman coeffient')

#Use whole matrix corrplot
# corrplot.mixed(Phi_cor_spearman, 
#                addCoef.col = 'black',
#                upper = "ellipse")
#Test for PerformanceAnalytics style correlationship plot
# library(PerformanceAnalytics)
# chart.Correlation(Phi_cor_spearman, histogram=TRUE, pch=19)

#Spearman PMI correlationship for by different cancer type-------------------------

#View(abs_word_pairs_PMI_sp)

# PMI_cor_spearman <- signif(cor(abs_word_pairs_PMI_sp[,-1],method = "spearman"),2)
# 
# corrplot(PMI_cor_spearman,
#          method = "color",
#          type = "upper",
#          addCoef.col = 'black',
#          title = ''
#          #title = 'PMI spearman coeffient by different cancer types'
#          )
# 
# text(1,1,'PMI spearman coeffient')
cor.m <- abs_word_pairs_PMI_sp[,-1]
colnames(cor.m) <- c("Bile duct" ,
                     "Colorectal" ,
                     "Esophageal" ,
                     "Liver" , 
                     "Pancreatic" ,
                     "Stomach")
PMI_cor_spearman <- signif(cor(cor.m,method = "spearman"),3)

corrplot(PMI_cor_spearman,
         method = "color",
         type = "upper",
         col = col(50),
         tl.pos = "d",
         mar=c(0,0,1.5,1),
         title = 'B. PMI spearman coeffient',
         cl.length = 11
)

corrplot(PMI_cor_spearman, 
         add = TRUE, 
         type = "lower", 
         method = "number", 
         diag = FALSE, 
         col = "black",
         tl.pos = "n", 
         mar=c(0,0,1.5,1),
         cl.pos = "n")
#Spearman Cosine correlationship for by different cancer type-------------------------

#View(abs_word_pairs_cosine_tfidf_sp)

# Cosine_cor_spearman <- signif(cor(abs_word_pairs_cosine_tfidf_sp[,-1],method = "spearman"),2)
# 
# corrplot(Cosine_cor_spearman,
#          method = "color",
#          type = "upper",
#          addCoef.col = 'black',
#          title = ''
#          #title = 'Cosine spearman coeffient by different cancer types'
#          )
# 
# text(1,1,'Cosine spearman coeffient')
cor.m <- abs_word_pairs_cosine_tfidf_sp[,-1]
colnames(cor.m) <- c("Bile duct" ,
                     "Colorectal" ,
                     "Esophageal" ,
                     "Liver" , 
                     "Pancreatic" ,
                     "Stomach")
Cosine_cor_spearman <- signif(cor(cor.m,method = "spearman"),3)

corrplot(Cosine_cor_spearman,
         method = "color",
         type = "upper",
         col = col(50),
         tl.pos = "d",
         mar=c(0,0,1.5,1),
         title = 'C. Cosine spearman coeffient',
         cl.length = 11
)

corrplot(Cosine_cor_spearman, 
         add = TRUE, 
         type = "lower", 
         method = "number", 
         diag = FALSE, 
         col = "black",
         tl.pos = "n", 
         mar=c(0,0,1.5,1),
         cl.pos = "n")
#Spearman Cosine correlationship for by different cancer type-------------------------

#View(abs_word_pairs_cosine_sp)

# Cosine_tfidf_cor_spearman <- signif(cor(abs_word_pairs_cosine_sp[,-1],method = "spearman"),2)
# 
# corrplot(Cosine_tfidf_cor_spearman,
#          method = "color",
#          type = "upper",
#          addCoef.col = 'black',
#          title = ''
#          #title = 'Cosine spearman coeffient by different cancer types'
# )
# 
# text(1,1,'Cosine tfidf spearman coeffient')
cor.m <- abs_word_pairs_cosine_sp[,-1]
colnames(cor.m) <- c("Bile duct" ,
                     "Colorectal" ,
                     "Esophageal" ,
                     "Liver" , 
                     "Pancreatic" ,
                     "Stomach")
Cosine_tfidf_cor_spearman <- signif(cor(cor.m,method = "spearman"),3)

corrplot(Cosine_tfidf_cor_spearman,
         method = "color",
         type = "upper",
         col = col(50),
         tl.pos = "d",
         mar=c(0,0,1.5,1),
         title = 'D. Cosine tfidf spearman coeffient',
         cl.length = 11
)

corrplot(Cosine_tfidf_cor_spearman, 
         add = TRUE, 
         type = "lower", 
         method = "number", 
         diag = FALSE, 
         col = "black",
         tl.pos = "n", 
         mar=c(0,0,1.5,1),
         cl.pos = "n")


#Change the par margin form of layout
par(mfrow = c(3,2))
#Spearman Bile duct cancer by three methods--------------------
#Use full_join to join
BTC_cor_spearman <- abs_word_pairs_phi_sp %>%
                    select(Gene,`Bile duct cancer`) %>%
                    full_join(abs_word_pairs_PMI_sp %>%
                                select(Gene,`Bile duct cancer`), by = 'Gene') %>%
                    full_join(abs_word_pairs_cosine_tfidf_sp %>%
                                select(Gene,`Bile duct cancer`), by = 'Gene') %>%
                    full_join(abs_word_pairs_cosine_sp %>%
                               select(Gene,`Bile duct cancer`), by = 'Gene') %>%
                    replace(is.na(.),0) 
colnames(BTC_cor_spearman) <- c('Genes','Phi','PMI','Cosine','Cosine_tfidf')

# BTC_cor_spearman <- signif(cor(BTC_cor_spearman[,-1],method = "spearman"),2)
# 
# corrplot(BTC_cor_spearman,
#          method = "color",
#          type = "upper",
#          addCoef.col = 'black',
#          title = ''
# )
# text(1,1,'Bile duct cancer spearman coeffient')
BTC_cor_spearman <- signif(cor(BTC_cor_spearman[,-1],method = "spearman"),3)

corrplot(BTC_cor_spearman,
         method = "color",
         type = "upper",
         col = col(50),
         tl.pos = "d",
         mar=c(0,0,1.5,1),
         title = 'A. Bile duct cancer spearman coeffient',
         cl.length = 11
)

corrplot(BTC_cor_spearman, 
         add = TRUE, 
         type = "lower", 
         method = "number", 
         diag = FALSE, 
         col = "black",
         tl.pos = "n", 
         mar=c(0,0,1.5,1),
         cl.pos = "n")

#Spearman Colorectal cancer by three methods--------------------
#Use full_join to join
COAD_cor_spearman <- abs_word_pairs_phi_sp %>%
  select(Gene,`Colorectal cancer`) %>%
  full_join(abs_word_pairs_PMI_sp %>%
              select(Gene,`Colorectal cancer`), by = 'Gene') %>%
  full_join(abs_word_pairs_cosine_tfidf_sp %>%
              select(Gene,`Colorectal cancer`), by = 'Gene') %>%
  full_join(abs_word_pairs_cosine_sp %>%
              select(Gene,`Colorectal cancer`), by = 'Gene') %>%
  replace(is.na(.),0) 
colnames(COAD_cor_spearman) <- c('Genes','Phi','PMI','Cosine','Cosine_tfidf')

# COAD_cor_spearman <- signif(cor(COAD_cor_spearman[,-1],method = "spearman"),2)
# 
# corrplot(COAD_cor_spearman,
#          method = "color",
#          type = "upper",
#          addCoef.col = 'black',
#          title = ''
# )
# text(1,1,'Colorectal cancer spearman coeffient')
COAD_cor_spearman <- signif(cor(COAD_cor_spearman[,-1],method = "spearman"),3)

corrplot(COAD_cor_spearman,
         method = "color",
         type = "upper",
         col = col(50),
         tl.pos = "d",
         mar=c(0,0,1.5,1),
         title = 'B. Colorectal cancer spearman coeffient',
         cl.length = 11
)

corrplot(COAD_cor_spearman, 
         add = TRUE, 
         type = "lower", 
         method = "number", 
         diag = FALSE, 
         col = "black",
         tl.pos = "n", 
         mar=c(0,0,1.5,1),
         cl.pos = "n")

#Spearman Esophageal cancer by three methods--------------------
#Use full_join to join
ESCA_cor_spearman <- abs_word_pairs_phi_sp %>%
  select(Gene,`Esophageal cancer`) %>%
  full_join(abs_word_pairs_PMI_sp %>%
              select(Gene,`Esophageal cancer`), by = 'Gene') %>%
  full_join(abs_word_pairs_cosine_tfidf_sp %>%
              select(Gene,`Esophageal cancer`), by = 'Gene') %>%
  full_join(abs_word_pairs_cosine_sp %>%
              select(Gene,`Esophageal cancer`), by = 'Gene') %>%
  replace(is.na(.),0) 
colnames(ESCA_cor_spearman) <- c('Genes','Phi','PMI','Cosine','Cosine_tfidf')

# ESCA_cor_spearman <- signif(cor(ESCA_cor_spearman[,-1],method = "spearman"),2)
# 
# corrplot(ESCA_cor_spearman,
#          method = "color",
#          type = "upper",
#          addCoef.col = 'black',
#          title = ''
# )
# text(1,1,'Esophageal cancer spearman coeffient')

ESCA_cor_spearman <- signif(cor(ESCA_cor_spearman[,-1],method = "spearman"),3)

corrplot(ESCA_cor_spearman,
         method = "color",
         type = "upper",
         col = col(50),
         tl.pos = "d",
         mar=c(0,0,1.5,1),
         title = 'C. Esophageal cancer spearman coeffient',
         cl.length = 11
)

corrplot(ESCA_cor_spearman, 
         add = TRUE, 
         type = "lower", 
         method = "number", 
         diag = FALSE, 
         col = "black",
         tl.pos = "n", 
         mar=c(0,0,1.5,1),
         cl.pos = "n")

#Spearman Liver cancer by three methods--------------------
#Use full_join to join
LIHC_cor_spearman <- abs_word_pairs_phi_sp %>%
  select(Gene,`Liver cancer`) %>%
  full_join(abs_word_pairs_PMI_sp %>%
              select(Gene,`Liver cancer`), by = 'Gene') %>%
  full_join(abs_word_pairs_cosine_tfidf_sp %>%
              select(Gene,`Liver cancer`), by = 'Gene') %>%
  full_join(abs_word_pairs_cosine_sp %>%
              select(Gene,`Liver cancer`), by = 'Gene') %>%
  replace(is.na(.),0) 
colnames(LIHC_cor_spearman) <- c('Genes','Phi','PMI','Cosine','Cosine_tfidf')

# LIHC_cor_spearman <- signif(cor(LIHC_cor_spearman[,-1],method = "spearman"),2)
# 
# corrplot(LIHC_cor_spearman,
#          method = "color",
#          type = "upper",
#          addCoef.col = 'black',
#          title = ''
# )
# text(1,1,'Liver cancer spearman coeffient')

LIHC_cor_spearman <- signif(cor(LIHC_cor_spearman[,-1],method = "spearman"),3)

corrplot(LIHC_cor_spearman,
         method = "color",
         type = "upper",
         col = col(50),
         tl.pos = "d",
         mar=c(0,0,1.5,1),
         title = 'D. Liver cancer spearman coeffient',
         cl.length = 11
)

corrplot(LIHC_cor_spearman, 
         add = TRUE, 
         type = "lower", 
         method = "number", 
         diag = FALSE, 
         col = "black",
         tl.pos = "n", 
         mar=c(0,0,1.5,1),
         cl.pos = "n")

#Spearman Pancreatic Cancer by three methods--------------------
#Use full_join to join
PAAD_cor_spearman <- abs_word_pairs_phi_sp %>%
  select(Gene,`Pancreatic Cancer`) %>%
  full_join(abs_word_pairs_PMI_sp %>%
              select(Gene,`Pancreatic Cancer`), by = 'Gene') %>%
  full_join(abs_word_pairs_cosine_tfidf_sp %>%
              select(Gene,`Pancreatic Cancer`), by = 'Gene') %>%
  full_join(abs_word_pairs_cosine_sp %>%
              select(Gene,`Pancreatic Cancer`), by = 'Gene') %>%
  replace(is.na(.),0) 
colnames(PAAD_cor_spearman) <- c('Genes','Phi','PMI','Cosine','Cosine_tfidf')

# PAAD_cor_spearman <- signif(cor(PAAD_cor_spearman[,-1],method = "spearman"),2)
# 
# corrplot(PAAD_cor_spearman,
#          method = "color",
#          type = "upper",
#          addCoef.col = 'black',
#          title = ''
# )
# text(1,1,'Pancreatic Cancer spearman coeffient')

PAAD_cor_spearman <- signif(cor(PAAD_cor_spearman[,-1],method = "spearman"),3)

corrplot(PAAD_cor_spearman,
         method = "color",
         type = "upper",
         col = col(50),
         tl.pos = "d",
         mar=c(0,0,1.5,1),
         title = 'E. Pancreatic Cancer spearman coeffient',
         cl.length = 11
)

corrplot(PAAD_cor_spearman, 
         add = TRUE, 
         type = "lower", 
         method = "number", 
         diag = FALSE, 
         col = "black",
         tl.pos = "n", 
         mar=c(0,0,1.5,1),
         cl.pos = "n")

#Spearman Stomach cancer by three methods--------------------
#Use full_join to join
STAD_cor_spearman <- abs_word_pairs_phi_sp %>%
  select(Gene,`Stomach cancer`) %>%
  full_join(abs_word_pairs_PMI_sp %>%
              select(Gene,`Stomach cancer`), by = 'Gene') %>%
  full_join(abs_word_pairs_cosine_tfidf_sp %>%
              select(Gene,`Stomach cancer`), by = 'Gene') %>%
  full_join(abs_word_pairs_cosine_sp %>%
              select(Gene,`Stomach cancer`), by = 'Gene') %>%
  replace(is.na(.),0) 
colnames(STAD_cor_spearman) <- c('Genes','Phi','PMI','Cosine','Cosine_tfidf')

# STAD_cor_spearman <- signif(cor(STAD_cor_spearman[,-1],method = "spearman"),2)
# 
# corrplot(STAD_cor_spearman,
#          method = "color",
#          type = "upper",
#          addCoef.col = 'black',
#          title = ''
# )
# text(1,1,'Stomach cancer spearman coeffient')

STAD_cor_spearman <- signif(cor(STAD_cor_spearman[,-1],method = "spearman"),3)

corrplot(STAD_cor_spearman,
         method = "color",
         type = "upper",
         col = col(50),
         tl.pos = "d",
         mar=c(0,0,1.5,1),
         title = 'F. Stomach cancer spearman coeffient',
         cl.length = 11
)

corrplot(STAD_cor_spearman, 
         add = TRUE, 
         type = "lower", 
         method = "number", 
         diag = FALSE, 
         col = "black",
         tl.pos = "n", 
         mar=c(0,0,1.5,1),
         cl.pos = "n")
# dev.off()

#Write to xlsx ------------------------------------------
library(xlsx)
filepath <- paste0('D:/wyq/消化道肿瘤热点基因文献挖掘/Table/','Spearman_coeffient.xlsx')
#?write.xlsx
filepath

write.xlsx(x = Phi_cor_spearman,
           file = filepath,
           sheetName = 'Phi_by_cancers',
           append = T)
write.xlsx(x = PMI_cor_spearman,
           file = filepath,
           sheetName = 'PMI_by_cancers',
           append = T)
write.xlsx(x = Cosine_cor_spearman,
           file = filepath,
           sheetName = 'Cosine_by_cancers',
           append = T)
write.xlsx(x = Cosine_tfidf_cor_spearman,
           file = filepath,
           sheetName = 'Cosine_tfidf_by_cancers',
           append = T)
write.xlsx(x = BTC_cor_spearman,
           file = filepath,
           sheetName = 'BTC_by_method',
           append = T)
write.xlsx(x = COAD_cor_spearman,
           file = filepath,
           sheetName = 'COAD_by_method',
           append = T)
write.xlsx(x = ESCA_cor_spearman,
           file = filepath,
           sheetName = 'ESCA_by_method',
           append = T)
write.xlsx(x = LIHC_cor_spearman,
           file = filepath,
           sheetName = 'LIHC_by_method',
           append = T)
write.xlsx(x = PAAD_cor_spearman,
           file = filepath,
           sheetName = 'PAAD_by_method',
           append = T)
write.xlsx(x = STAD_cor_spearman,
           file = filepath,
           sheetName = 'STAD_by_method',
           append = T)
