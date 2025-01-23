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
pro_data$replicate <- rep(c(1,2,3), times =14)
color <- c("#FDAE61","#ABDDA4","#FEE08B","#FE9929",
           "#C6DBEF","#9ECAE1","#6BAED6","#2171B5",
           "#9E9AC8","#6A51A3","#74C476","#006D2C","#8DD3C7","#F1B6DA")
pch <-c(16,20,16,16,15,15,15,15,18,23,24,17,20,20)
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
zlab<-paste0("PC2(",round(summ$importance[2,3]*100,2),"%)")
#形状：
#时期：15,16,17,18,20,21,22,23,24,25
x<- scatterplot3d(x = pro_data$PC1,y = pro_data$PC2,z = pro_data$PC3,
              xlab = 'PC1',ylab = 'PC2',zlab = 'PC3',cex.symbols = 3,
              color =color[factor(pro_data$tissue)],box = TRUE,grid = TRUE,
              pch =pch[factor(pro_data$tissue)],
              type = 'h',lty.grid = 2, 
              )
#addgrids3d(pro_data[, 1:3], grid = c("xy", "xz", "yz"),lty.grid =2,angle = 45 )
#add
zz.coords <- x$xyz.convert(pro_data[,1], pro_data[,2], pro_data[,3])


text(zz.coords$x,zz.coords$y,labels = pro_data$replicate,
   col="black",cex =0.6)

legend("topright", legend =c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                 "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2"),
           col = color, pch = pch,cex=0.4)
####新####






scatter3D(x = pro_data[,1], y = pro_data[,2], z =pro_data[,3], #bgvar = mag,
          pch = 21, cex = 3,col="black",bg=color,
          xlab = "PC1", ylab = "PC2",
          zlab = "PC3", 
          ticktype = "detailed",bty = "f",box = TRUE,
          #panel.first = panelfirst,
          theta = 45, phi = 45, d=3,
          colkey = FALSE) #list(length = 0.5, width = 0.5, cex.clab = 0.75))






####旧#####
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
rna_data$replicate <- rep(c(1,2,3), times =14)
color <- c("#FDAE61","#ABDDA4","#FEE08B","#FE9929",
           "#C6DBEF","#9ECAE1","#6BAED6","#2171B5",
           "#9E9AC8","#6A51A3","#74C476","#006D2C","#8DD3C7","#F1B6DA")
pch <-c(16,20,16,16,15,15,15,15,18,23,24,17,20,20)
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
zlab<-paste0("PC2(",round(summ$importance[2,3]*100,2),"%)")
#形状：
#时期：15,16,17,18,20,21,22,23,24,25
x<- scatterplot3d(x = rna_data$PC1,y = rna_data$PC2,z = rna_data$PC3,
                  xlab = 'PC1',ylab = 'PC2',zlab = 'PC3',cex.symbols = 3,
                  color =color[factor(rna_data$tissue)],box = TRUE,grid = TRUE,
                  pch =pch[factor(rna_data$tissue)],
                  type = 'h',lty.grid = 2, 
)
#addgrids3d(pro_data[, 1:3], grid = c("xy", "xz", "yz"),lty.grid =2,angle = 45 )
#add
zz.coords <- x$xyz.convert(rna_data[,1], rna_data[,2], rna_data[,3])


text(zz.coords$x,zz.coords$y,labels = rna_data$replicate,
     col="black",cex =0.6)

legend("topright", legend =c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                             "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2"),
       col = color, pch = pch,cex=0.4)

