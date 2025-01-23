####Spearman相关性与m6A level（平均水平）散点图+拟合线####
rm(list = ls())
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
####各个组织m6A基因蛋白与转录之间的Spearman相关性####
Tissue_info <-  c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")
C <- read.delim("../../../3_19/conserved.txt")
C <- unique(C[,13])
folder_path <-"D:/Workspace/Glycine_max/m6A/3_19/nonconserved/nonconserved/"
for (i in 1:14)
{
  # 构建文件路径
  file_path <- file.path(folder_path, paste0("non", i, "conserved.anno.peak.new.anno.txt"))
  # 导入数据
  Data <- read.delim(file_path)
  Data <- unique(Data[,14])
  name <- paste0("NC", i)
  assign(name, Data)
}
for (i in 1:14) {
  name <- paste0("NC", i)
  Data <- get(name)
  Data <- c(C,Data)
  Data <- unique(Data)
  name <- paste0("M", i)
  assign(name, Data)
}
##表达量关联##
PRO <- read.csv("D:/Workspace/Glycine_max/Data/PRO_Raw_Gene_Tissue.csv",row.names = 1)
###RNA <- read.csv("D:/Workspace/Glycine_max/Data/RNA_Tissue.csv",row.names = 1)###此处不要过滤后TPM文件
RNA <- read.csv("D:/Workspace/Glycine_max/Data/RNA_Tissue_1.csv",row.names = 1)
J <- intersect(rownames(PRO),rownames(RNA))
PRO <- subset(PRO,rownames(PRO) %in% J)
RNA <- subset(RNA,rownames(RNA) %in% J)
##Pairwise##
for (i in 1:nrow(PRO)) {
  for (j in 1:ncol(RNA)) {
    if (!is.na(PRO[i,j]) & is.na(RNA[i,j])) {PRO[i,j] <- NA}
    if (is.na(PRO[i,j]) & !is.na(RNA[i,j])) {RNA[i,j] <- NA}
  }
}
##准备结果文件##
Result_all <- data.frame(matrix("", nrow = 14, ncol =3))
colnames(Result_all) <- c("Tissue","m6A","non-m6A")
Result_all$Tissue <- Tissue_info

for (i in 1:14) 
{ 
  name <- paste0("M", i)
  data <-get(name)
  PRO_temp <- subset(PRO,rownames(PRO)%in%data)  ####第一个组织既有m6A又可能有蛋白表达水平的基因
  PRO_temp <- PRO_temp[!is.na(PRO_temp[,i]),]    ####第一个组织既有m6A又有蛋白表达水平的基因
  RNA_temp <-  subset(RNA,rownames(RNA)%in%data)####第一个组织既有m6A和蛋白表达水平，可能有的基因
  RNA_temp <- RNA_temp[!is.na(RNA_temp[,i]),]
  PRO_temp <-  subset(PRO_temp,rownames(PRO_temp)%in%rownames(RNA_temp))
  Result_all[i,2] <- cor(RNA_temp[,i],PRO_temp[,i],method = "spearman",use = "na.or.complete")
}
for (i in 1:14) 
{ 
  name <- paste0("M", i)
  data <-get(name)
  PRO_temp <- subset(PRO,!rownames(PRO)%in%data)
  PRO_temp <- PRO_temp[!is.na(PRO_temp[,i]),]
  RNA_temp <- subset(RNA,!rownames(RNA)%in%data)
  RNA_temp <- RNA_temp[!is.na(RNA_temp[,i]),]
  PRO_temp <-  subset(PRO_temp,rownames(PRO_temp)%in%rownames(RNA_temp))
  Result_all[i,3] <- cor(RNA_temp[,i],PRO_temp[,i],method = "spearman",use = "na.or.complete")
}
library(reshape2)

Data <- melt(Result_all,id.vars=c("Tissue"),
             measure.vars = c("m6A","non-m6A"),variable.name = "Spearman",value.name = "Cofficient")
Data$Cofficient <- as.numeric(Data$Cofficient)
Data$Tissue <- factor(Data$Tissue,levels = c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                             "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2"))
P_ALL <- ggplot(Data,aes(Spearman,Cofficient))+
  geom_point(aes(fill=Tissue),shape=21,size=4)+
  scale_fill_manual(values = c("#FDAE61","#ABDDA4","#FEE08B","#FE9929",
                               "#C6DBEF","#9ECAE1","#6BAED6","#2171B5",
                               "#9E9AC8","#6A51A3","#74C476","#006D2C","#8DD3C7","#F1B6DA"))+
  scale_y_continuous(expand = c(0,0),limits = c(0.15,0.6))+
  geom_line(aes(group=Tissue,color=Tissue),linewidth=0.6,linetype="dashed")+
  scale_color_manual(values = c("#FDAE61","#ABDDA4","#FEE08B","#FE9929",
                                "#C6DBEF","#9ECAE1","#6BAED6","#2171B5",
                                "#9E9AC8","#6A51A3","#74C476","#006D2C","#8DD3C7","#F1B6DA"))+
  theme_test()+ theme(legend.position = "none")+labs(x="")
####3'UTR####
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
####各个组织m6A基因蛋白与转录之间的Spearman相关性####
Tissue_info <-  c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")
C <- read.delim("../../../3_19/conserved.txt")
C <- C[C$annotation=="3' UTR",]
C <- unique(C[,13])
folder_path <-"D:/Workspace/Glycine_max/m6A/3_19/nonconserved/nonconserved/"
for (i in 1:14)
{
  # 构建文件路径
  file_path <- file.path(folder_path, paste0("non", i, "conserved.anno.peak.new.anno.txt"))
  # 导入数据
  Data <- read.delim(file_path)
  Data <-Data[Data$annotation=="3' UTR",]
  Data <- unique(Data[,14])
  name <- paste0("NC", i)
  assign(name, Data)
}
for (i in 1:14) {
  name <- paste0("NC", i)
  Data <- get(name)
  Data <- c(C,Data)
  Data <- unique(Data)
  name <- paste0("M", i)
  assign(name, Data)
}
##表达量关联##
PRO <- read.csv("D:/Workspace/Glycine_max/Data/PRO_VSN_Gene_Tissue.csv",row.names = 1)
###RNA <- read.csv("D:/Workspace/Glycine_max/Data/RNA_Tissue.csv",row.names = 1)###此处不要过滤后TPM文件
RNA <- read.csv("D:/Workspace/Glycine_max/Data/RNA_Tissue_1.csv",row.names = 1)




J <- intersect(rownames(PRO),rownames(RNA))
PRO <- subset(PRO,rownames(PRO) %in% J)
RNA <- subset(RNA,rownames(RNA) %in% J)
##Pairwise##
for (i in 1:nrow(PRO)) {
  for (j in 1:ncol(RNA)) {
    if (!is.na(PRO[i,j]) & is.na(RNA[i,j])) {PRO[i,j] <- NA}
    if (is.na(PRO[i,j]) & !is.na(RNA[i,j])) {RNA[i,j] <- NA}
  }
}
##准备结果文件##
Result_all <- data.frame(matrix("", nrow = 14, ncol =3))
colnames(Result_all) <- c("Tissue","m6A","non-m6A")
Result_all$Tissue <- Tissue_info

for (i in 1:14) 
{ 
  name <- paste0("M", i)
  data <-get(name)
  PRO_temp <- subset(PRO,rownames(PRO)%in%data)  ####第一个组织既有m6A又可能有蛋白表达水平的基因
  PRO_temp <- PRO_temp[!is.na(PRO_temp[,i]),]    ####第一个组织既有m6A又有蛋白表达水平的基因
  RNA_temp <-  subset(RNA,rownames(RNA)%in%data)####第一个组织既有m6A和蛋白表达水平，可能有的基因
  RNA_temp <- RNA_temp[!is.na(RNA_temp[,i]),]
  PRO_temp <-  subset(PRO_temp,rownames(PRO_temp)%in%rownames(RNA_temp))
  Result_all[i,2] <- cor(RNA_temp[,i],PRO_temp[,i],method = "spearman")
}
for (i in 1:14) 
{ 
  name <- paste0("M", i)
  data <-get(name)
  PRO_temp <- subset(PRO,!rownames(PRO)%in%data)
  PRO_temp <- PRO_temp[!is.na(PRO_temp[,i]),]
  RNA_temp <- subset(RNA,!rownames(RNA)%in%data)
  RNA_temp <- RNA_temp[!is.na(RNA_temp[,i]),]
  PRO_temp <-  subset(PRO_temp,rownames(PRO_temp)%in%rownames(RNA_temp))
  Result_all[i,3] <- cor(RNA_temp[,i],PRO_temp[,i],method = "spearman")
}
library(reshape2)

Data <- melt(Result_all,id.vars=c("Tissue"),
             measure.vars = c("m6A","non-m6A"),variable.name = "Spearman",value.name = "Cofficient")
Data$Cofficient <- as.numeric(Data$Cofficient)
Data$Tissue <- factor(Data$Tissue,levels = c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                             "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2"))
P_UTR <- ggplot(Data,aes(Spearman,Cofficient))+
  geom_point(aes(fill=Tissue),shape=21,size=4)+
  scale_fill_manual(values = c("#FDAE61","#ABDDA4","#FEE08B","#FE9929",
                               "#C6DBEF","#9ECAE1","#6BAED6","#2171B5",
                               "#9E9AC8","#6A51A3","#74C476","#006D2C","#8DD3C7","#F1B6DA"))+
  scale_y_continuous(expand = c(0,0),limits = c(0.15,0.6))+
  geom_line(aes(group=Tissue,color=Tissue),linewidth=0.6,linetype="dashed")+
  scale_color_manual(values = c("#FDAE61","#ABDDA4","#FEE08B","#FE9929",
                                "#C6DBEF","#9ECAE1","#6BAED6","#2171B5",
                                "#9E9AC8","#6A51A3","#74C476","#006D2C","#8DD3C7","#F1B6DA"))+
  theme_test()+ theme(legend.position = "none")+
  labs(y="")
####CDS####
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
####各个组织m6A基因蛋白与转录之间的Spearman相关性####
Tissue_info <-  c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")
C <- read.delim("../../../3_19/conserved.txt")
C <- C[grep("exon", C$annotation), ]
C <- unique(C[,13])
folder_path <-"D:/Workspace/Glycine_max/m6A/3_19/nonconserved/nonconserved/"
for (i in 1:14)
{
  # 构建文件路径
  file_path <- file.path(folder_path, paste0("non", i, "conserved.anno.peak.new.anno.txt"))
  # 导入数据
  Data <- read.delim(file_path)
  Data <-Data[grep("exon", Data$annotation), ]
  Data <- unique(Data[,14])
  name <- paste0("NC", i)
  assign(name, Data)
}
for (i in 1:14) {
  name <- paste0("NC", i)
  Data <- get(name)
  Data <- c(C,Data)
  Data <- unique(Data)
  name <- paste0("M", i)
  assign(name, Data)
}
##表达量关联##
PRO <- read.csv("D:/Workspace/Glycine_max/Data/PRO_VSN_Gene_Tissue.csv",row.names = 1)
###RNA <- read.csv("D:/Workspace/Glycine_max/Data/RNA_Tissue.csv",row.names = 1)###此处不要过滤后TPM文件
RNA <- read.csv("D:/Workspace/Glycine_max/Data/RNA_Tissue_1.csv",row.names = 1)




J <- intersect(rownames(PRO),rownames(RNA))
PRO <- subset(PRO,rownames(PRO) %in% J)
RNA <- subset(RNA,rownames(RNA) %in% J)
##Pairwise##
for (i in 1:nrow(PRO)) {
  for (j in 1:ncol(RNA)) {
    if (!is.na(PRO[i,j]) & is.na(RNA[i,j])) {PRO[i,j] <- NA}
    if (is.na(PRO[i,j]) & !is.na(RNA[i,j])) {RNA[i,j] <- NA}
  }
}
##准备结果文件##
Result_all <- data.frame(matrix("", nrow = 14, ncol =3))
colnames(Result_all) <- c("Tissue","m6A","non-m6A")
Result_all$Tissue <- Tissue_info

for (i in 1:14) 
{ 
  name <- paste0("M", i)
  data <-get(name)
  PRO_temp <- subset(PRO,rownames(PRO)%in%data)  ####第一个组织既有m6A又可能有蛋白表达水平的基因
  PRO_temp <- PRO_temp[!is.na(PRO_temp[,i]),]    ####第一个组织既有m6A又有蛋白表达水平的基因
  RNA_temp <-  subset(RNA,rownames(RNA)%in%data)####第一个组织既有m6A和蛋白表达水平，可能有的基因
  RNA_temp <- RNA_temp[!is.na(RNA_temp[,i]),]
  PRO_temp <-  subset(PRO_temp,rownames(PRO_temp)%in%rownames(RNA_temp))
  Result_all[i,2] <- cor(RNA_temp[,i],PRO_temp[,i],method = "spearman")
}
for (i in 1:14) 
{ 
  name <- paste0("M", i)
  data <-get(name)
  PRO_temp <- subset(PRO,!rownames(PRO)%in%data)
  PRO_temp <- PRO_temp[!is.na(PRO_temp[,i]),]
  RNA_temp <- subset(RNA,!rownames(RNA)%in%data)
  RNA_temp <- RNA_temp[!is.na(RNA_temp[,i]),]
  PRO_temp <-  subset(PRO_temp,rownames(PRO_temp)%in%rownames(RNA_temp))
  Result_all[i,3] <- cor(RNA_temp[,i],PRO_temp[,i],method = "spearman")
}
library(reshape2)

Data <- melt(Result_all,id.vars=c("Tissue"),
             measure.vars = c("m6A","non-m6A"),variable.name = "Spearman",value.name = "Cofficient")
Data$Cofficient <- as.numeric(Data$Cofficient)
Data$Tissue <- factor(Data$Tissue,levels = c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                             "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2"))
P_CDS<- ggplot(Data,aes(Spearman,Cofficient))+
  geom_point(aes(fill=Tissue),shape=21,size=4)+
  scale_fill_manual(values = c("#FDAE61","#ABDDA4","#FEE08B","#FE9929",
                               "#C6DBEF","#9ECAE1","#6BAED6","#2171B5",
                               "#9E9AC8","#6A51A3","#74C476","#006D2C","#8DD3C7","#F1B6DA"))+
  scale_y_continuous(expand = c(0,0),limits = c(0.2,0.65))+
  geom_line(aes(group=Tissue,color=Tissue),linewidth=0.6,linetype="dashed")+
  scale_color_manual(values = c("#FDAE61","#ABDDA4","#FEE08B","#FE9929",
                                "#C6DBEF","#9ECAE1","#6BAED6","#2171B5",
                                "#9E9AC8","#6A51A3","#74C476","#006D2C","#8DD3C7","#F1B6DA"))+
  theme_test()+
  labs(y="")+labs(x="")
  
library(patchwork)
P_ALL+P_UTR+P_CDS+plot_layout(nrow=1)



####散点图####

