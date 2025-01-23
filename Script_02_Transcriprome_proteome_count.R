####转录组与蛋白组渐变色柱子####
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(paletteer)
library(ggpubr)
library(scales)
library(patchwork)
library(reshape2)
library(ggthemes)
##color##
RNA_COL=paletteer_c("ggthemes::Orange-Blue-White Diverging", n =14)
PRO_COL=paletteer_c("ggthemes::Red-Green-White Diverging", n=14)
##load data##
PRO <-read.csv("D:/Workspace/Glycine_max/Data/PRO_Raw_Gene_Tissue.csv",row.names = 1)
RNA <-read.csv("D:/Workspace/Glycine_max/Data/RNA_Tissue.csv",row.names = 1)
PRO$Count <- apply(PRO, 1, function(row) sum(!is.na(row)))
RNA$Count <- apply(RNA, 1, function(row) sum(!is.na(row)))
for (i in 1:nrow(PRO)) {
  non_na_indices <- !is.na(PRO[i, 1:14])
  PRO[i, 1:14][non_na_indices] <- PRO[i, 15]
}
for (i in 1:nrow(RNA)) {
  non_na_indices <- !is.na(RNA[i, 1:14])
  RNA[i, 1:14][non_na_indices] <- RNA[i, 15]
}
#Convert data#
PRO_Data <- melt(PRO[,-15])
RNA_Data <- melt(RNA[,-15])
colnames(PRO_Data)<- c("Tissue","Frequency")
colnames(RNA_Data)<- c("Tissue","Frequency")
PRO_Data <- PRO_Data[!(is.na(PRO_Data[, 2])), ]
RNA_Data <- RNA_Data[!(is.na(RNA_Data[, 2])), ]
RNA_Data$Frequency <- factor(RNA_Data$Frequency)
RNA_Data$Tissue <- factor(RNA_Data$Tissue,levels =rev(c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                                            "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2")))
RNA1 <- ggplot(RNA_Data,aes(x=Tissue))+
  geom_bar(aes(fill=Frequency),width = 0.7)+
  geom_hline(aes(yintercept=sum(RNA[, 15] == 14)), colour="black", linetype = "dashed", linewidth= 1)+
  scale_fill_manual(values =RNA_COL)+
  theme_half_open()+
  coord_flip()+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_text(size = 8))+
  theme(legend.text = element_text(size = 8))+
  theme(legend.key.size = unit(0.1, "cm"))  + theme(legend.title = element_text(size = 8))+
  labs(x = NULL, y = NULL)+
  theme(legend.position = "none")
RNA1 

PRO_Data$Frequency <- factor(PRO_Data$Frequency)
PRO_Data$Tissue <- factor(PRO_Data$Tissue,levels =rev(c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                                            "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2")))
PRO1 <- ggplot(PRO_Data,aes(x=Tissue))+
  geom_bar(aes(fill=Frequency),width = 0.7)+
  geom_hline(aes(yintercept=sum(PRO[, 15] == 14)), colour="black", linetype = "dashed", linewidth= 1)+
  scale_fill_manual(values =PRO_COL)+
  theme_half_open()+
  coord_flip()+
  ylim(0,15000)+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.text.y = element_text(size = 8))+
  theme(legend.text = element_text(size = 8))+
  theme(legend.key.size = unit(0.1, "cm"))  + theme(legend.title = element_text(size = 8))+
  labs(x = NULL, y = NULL)+
  theme(legend.position = "none")
PRO1
PRO1 + RNA1 + plot_layout(guides = "collect", widths = c(1, 1))
