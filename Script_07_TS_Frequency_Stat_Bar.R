###Frequency + Stat + Category
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(cowplot)
Frequency_adjust <- data.frame(matrix(NA, nrow =12758, ncol = 14))
PRO_co_quantified_sample <- read.csv("./Hou/PRO/PRO_co_quantified_sample.csv",row.names = 1)
PRO_Score_Tissue <- read.csv(file = "./Hou/PRO/PRO_Score_tissue.csv",row.names = 1)
PRO_Specific <- read.csv("./Hou/PRO/PRO_Specific.csv",row.names = 1)
PRO_Enriched_not_specific <-  read.csv("./Hou/PRO/PRO_Enriched_not_specific.csv",row.names = 1)
PRO_housekeeping <-  read.csv("./Hou/PRO/PRO_housekeeping.csv",row.names = 1)
PRO_others <- read.csv("./Hou/PRO/PRO_others.csv",row.names = 1)
RNA_Specific <- read.csv("./Hou/RNA/RNA_Specific.csv",row.names = 1)
RNA_Enriched_not_specific <-  read.csv("./Hou/RNA/RNA_Enriched_not_specific.csv",row.names = 1)
RNA_housekeeping <-  read.csv("./Hou/RNA/RNA_housekeeping.csv",row.names = 1)
RNA_others <- read.csv("./Hou/RNA/RNA_others.csv",row.names = 1)
for (i in 1:nrow(PRO_co_quantified_sample)) {
  # adjust NA
  for (j in 1:14) {
    cols <- ((j - 1) * 3 + 1):((j - 1) * 3 + 3)
    if (any(!is.na(PRO_co_quantified_sample[i, cols]))) {
      Frequency_adjust[i, j] <- 1
    }
  }
}
rownames(Frequency_adjust) <- rownames(PRO_Score_Tissue)
colnames(Frequency_adjust) <- colnames(PRO_Score_Tissue)
Frequency_adjust$Frequency <- rowSums(Frequency_adjust, na.rm = TRUE)
Data <- data.frame(matrix(NA, nrow =12758, ncol = 2))
colnames(Data) <- c("Frequency","Category")
rownames(Data) <- rownames(Frequency_adjust)
Data$Frequency <- Frequency_adjust$Frequency
for (i in 1:nrow(Data)) {
  
  if (rownames(Data)[i] %in% rownames(PRO_Specific)) {
    
    Data[i, 2] <- "Tissue Specific"
  }
}
for (i in 1:nrow(Data)) {
  
  if (rownames(Data)[i] %in% rownames(PRO_Enriched_not_specific)) {
 
    Data[i, 2] <- "Tissue Enriched not Specific"
  }
}
for (i in 1:nrow(Data)) {

  if (rownames(Data)[i] %in% rownames(PRO_others)) {

    Data[i, 2] <- "Others"
  }
}
for (i in 1:nrow(Data)) {

  if (rownames(Data)[i] %in% rownames(PRO_housekeeping)) {

    Data[i, 2] <- "HouseKeeping"
  }
}
Data$Category <- factor(Data$Category,levels = c("Tissue Specific","Tissue Enriched not Specific","HouseKeeping","Others"))
P1 <- ggplot(data=Data,aes(x=Frequency))+
  geom_bar(aes(fill=Category),
           position = 'stack',
           color='black',width=0.8,linewidth=0.4)+
  labs(y=NULL)+
  scale_fill_manual(values = c("#E41A1C","#BEBADA","#6BAED6","#8DD3C7"))+
  theme_half_open()+
  theme(legend.position = c(0.1,0.9))+
  theme(legend.key.size = unit(10, "pt"))+ 
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(colour = NA))+ 

  theme(legend.background = element_rect(fill = NA))+
  scale_x_continuous(breaks=seq(1,14,1), labels=seq(1,14)) +
  scale_y_continuous(breaks=seq(0,7000,1000), labels=c(0,1000,2000,3000,4000,5000,6000,7000))+
  labs(x="number of observed tissues")

####plot####
#protein count data#
PRO_Stat_All <- data.frame(matrix("", nrow = 4, ncol =3))
colnames(PRO_Stat_All) <- c("Type","Category","Number")
PRO_Stat_All$Type <- "Proteome"
PRO_Stat_All[1,2] <- "Tissue Eniched not Specific"
PRO_Stat_All[1,3] <- nrow(PRO_Enriched_not_specific)
PRO_Stat_All[2,2] <- "Tissue Specific"
PRO_Stat_All[2,3] <- nrow(PRO_Specific)
PRO_Stat_All[3,2] <- "HouseKeeping"
PRO_Stat_All[3,3] <- nrow(PRO_housekeeping)
PRO_Stat_All[4,2] <- "Others"
PRO_Stat_All[4,3] <- nrow(PRO_others)
#RNA count data
RNA_Stat_All <- data.frame(matrix("", nrow = 4, ncol =3))
colnames(RNA_Stat_All) <- c("Type","Category","Number")
RNA_Stat_All$Type <- "Transcriptome"
RNA_Stat_All[1,2] <- "Tissue Eniched not Specific"
RNA_Stat_All[1,3] <- nrow(RNA_Enriched_not_specific)
RNA_Stat_All[2,2] <- "Tissue Specific"
RNA_Stat_All[2,3] <- nrow(RNA_Specific)
RNA_Stat_All[3,2] <- "HouseKeeping"
RNA_Stat_All[3,3] <- nrow(RNA_housekeeping)
RNA_Stat_All[4,2] <- "Others"
RNA_Stat_All[4,3] <- nrow(RNA_others)
Stat <- rbind(PRO_Stat_All,RNA_Stat_All)
Stat$Number <- as.numeric(Stat$Number)
Stat$Type<- factor(Stat$Type,levels=c("Transcriptome","Proteome"))
Stat$Category <- factor(Stat$Category,levels =rev(c("Tissue Specific","Tissue Eniched not Specific","HouseKeeping","Others")))
P2 <- ggplot(data=Stat)+
  geom_col(aes(x=Type,y=Number,fill=Category),width=.5)+
  theme_classic()+
  scale_fill_manual(values = c("#8DD3C7","#6BAED6","#BEBADA","#E41A1C"))+
  theme(legend.key.size = unit(10, "pt"))+ 
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(colour = NA))+ 
  theme(legend.background = element_rect(fill = NA))+
  theme(legend.position ="none") +
  scale_y_continuous(breaks=seq(0,10000,2500), labels=c(0,2500,5000,7500,10000))+
  labs(x=NULL)+
coord_flip()
library(patchwork)
P1/P2+plot_layout(heights =  c(3, 1))
