####准备软件包####
rm(list = ls())
library(ggplot2)
library(plyr)
library(scales)
library(ggcorrplot)
library(grid)
library(ggbiplot)
library(tidyverse)
library(ggstar)
library(paletteer)
library(ggrepel)
library(RColorBrewer)
library(factoextra)
library(FactoMineR)
library(cowplot)
library(scatterplot3d)
library(ComplexHeatmap)
library(FactoClass)
####PRO水平####
gene_expression <- read.csv("D:/Workspace/Glycine_max/Data/PRO_VSN_Gene_Sample.csv", row.names=1)
Sample_info_sample <- read.csv("D:/Workspace/Glycine_max/Tissue_Specifity/Tissue_Specifity/Hou/Sample_info_Sample.csv")
gene_expression=as.matrix(gene_expression)
gene_expression_t=t(gene_expression)
gene_expression_filter=gene_expression_t[,which(colSums(gene_expression_t)>0)]
dim(gene_expression_filter)
w.pca<- prcomp(gene_expression_filter,scale. = TRUE)
df1 <- data.frame(w.pca$x, Sample_info_sample$Tissue)
colnames(df1)[43] <- "Tissue"
summ<-summary(w.pca)
pro_data <- as.data.frame(w.pca$x)
pro_data <- pro_data[,1:3]
pro_data$sample <- rownames(pro_data)
pro_data$tissue <- Sample_info_sample$Tissue
pro_data$tissue <- factor(pro_data$tissue,levels =c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                                         "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2"))
##因子排序
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
##ggplot2画图
##
df1$Tissue <- factor(df1$Tissue,levels =c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                          "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2"))

ggplot(data = df1,aes(x=PC1,y=PC2))+
  geom_star(size=4,aes(fill=Tissue,starshape=Tissue),color="black")+
  scale_fill_manual(values = c("#FDAE61","#ABDDA4","#FEE08B","#FE9929",
                               "#C6DBEF","#9ECAE1","#6BAED6","#2171B5",
                               "#9E9AC8","#6A51A3","#74C476","#006D2C","#8DD3C7","#F1B6DA"))+
  scale_starshape_manual(values =c(25,15,15,15,11,11,11,11,12,12,14,14,5,6))+
  #  geom_text_repel(aes(PC1, PC2,label=rownames(df1)),size =3 )+
  labs(x=xlab,y =ylab)+
  geom_hline(aes(yintercept=0), colour="black", linetype="dashed")+
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed")+
  theme_test()

####RNA水平####
rm(list = ls())
gene_expression <- read.csv("D:/Workspace/Glycine_max/Data/RNA_Sample_1.csv", row.names=1)
#样本信息#
Sample_info_sample <- read.csv("D:/Workspace/Glycine_max/Tissue_Specifity/Tissue_Specifity/Hou/Sample_info_Sample.csv")
gene_expression=as.matrix(gene_expression)
gene_expression_t=t(gene_expression)
gene_expression_filter=gene_expression_t[,which(colSums(gene_expression_t)>0)]
dim(gene_expression_filter)
w.pca<- prcomp(gene_expression_filter,scale. = TRUE)
df1 <- data.frame(w.pca$x, Sample_info_sample$Tissue)
colnames(df1)[43] <- "Tissue"
summ<-summary(w.pca)
rna_data <- as.data.frame(w.pca$x)
rna_data <- rna_data[,1:3]
rna_data$sample <- rownames(rna_data)
rna_data$tissue <- Sample_info_sample$Tissue
rna_data$tissue <- factor(rna_data$tissue,levels =c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                                    "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2"))
##因子排序
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
##ggplot2画图
##
df1$Tissue <- factor(df1$Tissue,levels =c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                          "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2"))

ggplot(data = df1,aes(x=PC1,y=PC2))+
  geom_star(size=4,aes(fill=Tissue,starshape=Tissue),color="black")+
  scale_fill_manual(values = c("#FDAE61","#ABDDA4","#FEE08B","#FE9929",
                               "#C6DBEF","#9ECAE1","#6BAED6","#2171B5",
                               "#9E9AC8","#6A51A3","#74C476","#006D2C","#8DD3C7","#F1B6DA"))+
  scale_starshape_manual(values =c(25,15,15,15,11,11,11,11,12,12,14,14,5,6))+
  #  geom_text_repel(aes(PC1, PC2,label=rownames(df1)),size =3 )+
  labs(x=xlab,y =ylab)+
  geom_hline(aes(yintercept=0), colour="black", linetype="dashed")+
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed")+
  theme_test()
