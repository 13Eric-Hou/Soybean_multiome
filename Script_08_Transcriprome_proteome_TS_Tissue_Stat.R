####各个组织RNA Speceifc & RNA Enriched not Specific数量
###PRO Speceifc & PRO Enriched not Specific
##准备文件：PRO_Enriched_not_specfic    PRO_Specific
##          RNA_Enriched_not_specfic    RNA_Specific
rm(list = ls())
PRO_Enriched_not_specific <- read.csv("D:/Workspace/Glycine_max/Tissue_Specifity/Tissue_Specifity/Hou/PRO/PRO_Enriched_not_specific.csv",row.names = 1)
PRO_Specific <- read.csv("D:/Workspace/Glycine_max/Tissue_Specifity/Tissue_Specifity/Hou/PRO/PRO_Specific.csv",row.names = 1)
RNA_Enriched_not_specific <- read.csv("D:/Workspace/Glycine_max/Tissue_Specifity/Tissue_Specifity/Hou/RNA/RNA_Enriched_not_specific.csv",row.names = 1)
RNA_Specific <- read.csv("D:/Workspace/Glycine_max/Tissue_Specifity/Tissue_Specifity/Hou/RNA/RNA_Specific.csv",row.names = 1)
PRO_Score_Tissue<-read.csv("D:/Workspace/Glycine_max/Tissue_Specifity/Tissue_Specifity/Hou/PRO/PRO_Score_tissue.csv",row.names = 1)
RNA_Score_Tissue<-read.csv("D:/Workspace/Glycine_max/Tissue_Specifity/Tissue_Specifity/Hou/RNA/RNA_Score_tissue.csv",row.names = 1)

PRO_Specific_Stat <- data.frame(matrix("", nrow = 14, ncol =4))
colnames(PRO_Specific_Stat) <- c("Tissue","Number","Type","Category")
PRO_Specific_Stat$Tissue <- colnames(PRO_Score_Tissue)
PRO_Specific_Stat$Type <- "PRO"
PRO_Specific_Stat$Category <- "PRO Specific"
PRO_Enriched_not_specific_Stat <- data.frame(matrix("", nrow = 14, ncol =4))
colnames(PRO_Enriched_not_specific_Stat) <- c("Tissue","Number","Type","Category")
PRO_Enriched_not_specific_Stat$Tissue <- colnames(PRO_Score_Tissue)
PRO_Enriched_not_specific_Stat$Type <- "PRO"
PRO_Enriched_not_specific_Stat$Category <- "PRO Enriched not Specific"
RNA_Specific_Stat <- data.frame(matrix("", nrow = 14, ncol =4))
colnames(RNA_Specific_Stat) <- c("Tissue","Number","Type","Category")
RNA_Specific_Stat$Tissue <- colnames(PRO_Score_Tissue)
RNA_Specific_Stat$Type <- "RNA"
RNA_Specific_Stat$Category <- "RNA Specific"
RNA_Enriched_not_specific_Stat <- data.frame(matrix("", nrow = 14, ncol =4))
colnames(RNA_Enriched_not_specific_Stat) <- c("Tissue","Number","Type","Category")
RNA_Enriched_not_specific_Stat$Tissue <- colnames(PRO_Score_Tissue)
RNA_Enriched_not_specific_Stat$Type <- "RNA"
RNA_Enriched_not_specific_Stat$Category <- "RNA Enriched not Specific"
for (i in 1:14) {
  # 统计"PRO_enriched"单元格数量
  PRO_Specific_Stat[i, 2] <- sum(PRO_Specific[, i] >= 4, na.rm = TRUE)
  PRO_Enriched_not_specific_Stat[i, 2] <- sum(PRO_Enriched_not_specific[, i] >=2.5, na.rm = TRUE)
  RNA_Specific_Stat[i, 2] <- sum(RNA_Specific[, i] >= 4, na.rm = TRUE)
  RNA_Enriched_not_specific_Stat[i, 2] <- sum(RNA_Enriched_not_specific[, i] >= 2.5, na.rm = TRUE)
}
Stat_Data <- rbind(PRO_Enriched_not_specific_Stat,PRO_Specific_Stat,RNA_Enriched_not_specific_Stat,RNA_Specific_Stat)
Stat_Data$Number <- as.numeric(Stat_Data$Number)
Stat_Data$Category <-factor(Stat_Data$Category,levels=c("PRO Specific","PRO Enriched not Specific",
                                                        "RNA Specific", "RNA Enriched not Specific"))

Stat_Data$Tissue= factor(Stat_Data$Tissue,levels =rev(c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                                    "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2")))
##绘图
library(RColorBrewer)
library(tidyverse)
library(cowplot)
library(ggthemes)
ggplot(Stat_Data, aes(x=Tissue)) +
  geom_bar(data=subset(Stat_Data,Type=="PRO"), aes(y=Number, fill=Category), stat="identity") +
  geom_bar(data=subset(Stat_Data,Type=="RNA"), aes(y=-Number, fill=Category), stat="identity") +
  geom_hline(yintercept=0, colour="white", lwd=2) +
  theme_half_open()+
  coord_flip(ylim=c(-1000,1000)) + 
  scale_fill_manual(values = c("#A1D99B","#238B45","#9ECAE1","#2171B5"))+
  scale_y_continuous(breaks=seq(-1000,1000,250), labels=c(1000,750,500,250,0,250,500,750,1000)) +
  labs(y="Number", x="Tissue") +
  ggtitle("RNA | PRO") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.02, 0.99),legend.key.size = unit(7, "pt"))+ 
  theme(legend.text = element_text(size = 7),
        legend.title = element_text(colour = NA))+ 
  theme(legend.background = element_rect(fill = NA))
