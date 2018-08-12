#Follow the instructions below
#https://github.com/js229/Vennerable
#source("http://bioconductor.org/biocLite.R")
#install.packages(c("graph", "RBGL", "xtable"))
#biocLite(c("graph", "RBGL"))
library(devtools)
#install_github("js229/Vennerable")
library(tidyverse)
library(Vennerable)
library(Cairo)
#First step to get the PMID of 6 GI cancer literauture

GiPath <- 'D:/wyq/Gi_Gene_Timeline/'

setwd(GiPath)

#Get PMID by list.files
CholPmid <- list.files('./chol_abstract/') %>% str_remove_all('.txt')
#Test for correct
#CholPmid <-str_remove_all(CholPmid,'.txt')
#test <- str_remove_all(CholPmid,'.txt')
CocaPmid <- list.files('./coca_abstract/') %>% str_remove_all('.txt')
EscaPmid <- list.files('./esca_abstract/') %>% str_remove_all('.txt')
LihcPmid <- list.files('./lihc_abstract/') %>% str_remove_all('.txt')
PaadPmid <- list.files('./paad_abstract/') %>% str_remove_all('.txt')
StadPmid <- list.files('./stad_abstract/') %>% str_remove_all('.txt')


#Plot vennerable 
data<-Venn(list("胆管癌"=CholPmid,
                "肠癌"=CocaPmid,
                "食道癌"=EscaPmid,
                "肝癌"=LihcPmid,
                "胰腺癌"=PaadPmid,
                "胃癌"=StadPmid))
cairo_pdf('Venn_6_GI_cancer.pdf',family = 'GB1')
plot(data,doWeight=F)
# #Add legends
# legend("bottomright", 
#        legend = c("Group 1", "Group 2")#, 
#        # col = c(rgb(0.2,0.4,0.1,0.7), 
#        #         rgb(0.8,0.4,0.1,0.7)), 
#        # pch = c(17,19), 
#        # bty = "n", 
#        # pt.cex = 2, 
#        # cex = 1.2, 
#        # text.col = "black", 
#        # horiz = F , 
#        # inset = c(0.1, 0.1)
#        )
dev.off()

VennTable<-data@IndicatorWeight
#str(VennTable)
VennTable <- as.tibble(VennTable) %>% arrange(desc(`.Weight`))
write.csv(VennTable,file = 'VennTableInfo.csv')


                                   