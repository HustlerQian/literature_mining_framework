#Use for annotation for precision and recall
#Load template
#library(xlsx)
# chol=unique(chol_gene.df$Concept)
# coca=unique(coca_gene.df$Concept)
# esca=unique(esca_gene.df$Concept)
# lihc=unique(lihc_gene.df$Concept)
# 
# paad_gene.df <- read.delim('./paad_metamap/conceptsInFile_Gene.txt',sep = '\t')
# 
# colnames(paad_gene.df) <- c('PMID','Date','Concept','Count')
# paad=unique(paad_gene.df$Concept)
# stad=unique(stad_gene.df$Concept)


gene_list=unique(GI.tdf$Concept)

#gongci
gongci=data.frame(cooccurence=gene_list)
write.csv(gongci,file = 'gongci.csv',row.names = F)


#Find regular error
Regular.tdf <- bind_rows(stad_gene2.df,paad_gene2.df,lihc_gene2.df,esca_gene2.df,coca_gene2.df,chol_gene2.df)


regular_gene_list <- unique(Regular.tdf$Concept)

regular<-data.frame(regular=regular_gene_list)

write.csv(regular,file = 'regular.csv',row.names = F)
