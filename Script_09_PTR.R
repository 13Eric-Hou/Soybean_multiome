####PTR计算####
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggcor)
library(ggthemes)
##需要的文件：蛋白表达量；TPM表达量；（TPM小于1的表达量视为NA）
##蛋白表达量/TPM表达量的结果，再取以2为底的对数；即是PTR
##数据文件导入##
rm(list = ls())
PRO <- read.csv("D:/Workspace/Glycine_max/Data/PRO_Raw_Gene_Tissue.csv",row.names = 1)
RNA <- read.csv("D:/Workspace/Glycine_max/Data/RNA_Tissue_1.csv",row.names = 1)
PRO <- subset(PRO,rownames(PRO)%in%rownames(RNA))
RNA <- subset(RNA,rownames(RNA)%in%rownames(PRO))
##循环保留同时有的##
PRO[!is.na(PRO) & is.na(RNA)]<- NA
RNA[!is.na(RNA) & is.na(PRO)]<- NA
identical(is.na(RNA), is.na(PRO))    ##判断是否转换成功
##计算PTR##
PTR <-data.frame(matrix(nrow = 12855, ncol = 14))
colnames(PTR) <-colnames(PRO)
rownames(PTR) <-rownames(RNA)
for (i in 1:nrow(PRO))
{
  for (j in 1:ncol(PRO))
  {
    if (!is.na(PRO[i,j])&!is.na(RNA[i,j])) {
      PTR[i,j] <- log2(PRO[i,j]/RNA[i,j])
      
    }
  }
}

write.csv(PTR,file = "./PTR_1.csv")
##筛选出至少具有10对成对观测值的基因##
PRO_n<-PRO[rowSums(is.na(PRO)) <= 4, ]
RNA_n <-RNA[rowSums(is.na(RNA)) <= 4, ]
PTR_n <- subset(PTR,rownames(PTR) %in% rownames(PRO_n))
##准备曲线图文件
Data <- data.frame(matrix(nrow = 9008, ncol = 15))
colnames(Data) <- c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8","Median")
rownames(Data) <- rownames(PTR_n)
for (i in 1:14) {
  Data[,i] <- PTR_n[,i]
}
##Median PTR统计##
Data[,15] <-apply(Data[, 1:14], 1, median, na.rm = TRUE)
median <- median(Data$Median)
sd <- sd(Data$Median)
##画图##
#整体的图
data <- melt(Data[,c(1:14)])
data <- data[!is.na(data[, 2]), ]
data$variable <- factor(data$variable,levels =c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                                "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2"))
ggplot(data,aes(value))+
  geom_density(aes(color=variable),linewidth=0.6)+
  geom_vline(xintercept = median+2*sd,linetype = "dashed", linewidth= 0.5)+
  geom_vline(xintercept = median-2*sd,linetype = "dashed", linewidth = 0.5)+
  scale_color_manual(values = c("#FDAE61","#ABDDA4","#FEE08B","#FE9929",
                                    "#C6DBEF","#9ECAE1","#6BAED6","#2171B5",
                                    "#9E9AC8","#6A51A3","#74C476","#006D2C","#8DD3C7","#F1B6DA"))+
  scale_x_continuous(breaks=seq(-10,10,5), labels=c(-10,-5,0,5,10))+
  theme_half_open()+
  theme(legend.key.size = unit(8,"pt"))+
 theme(legend.text = element_text(size = 8),
    legend.position = c(0.1, 0.76)) + theme(legend.title = element_text(colour = NA))
#median PTR#
data <- melt(Data[,15])
dense = data.frame(density(data$value,adjust = 1)[c('x','y')])
ggplot(data,aes(value))+
  geom_density(adjust=1,linewidth=1,fill="#FDAE6B",color="#FDAE6B")+
  geom_area(data = subset(dense,x<=(median-2*sd)), aes(x, y),fill="#B288B6",color="#B288B6")+
  geom_area(data = subset(dense,x>=(median+2*sd)), aes(x, y),fill="#41ABE2",color="#41ABE2")+
  geom_boxplot(color="#FDAE6B",width=0.007,position = position_nudge(x=0, y=-0.005), outlier.size = 1)+
  geom_vline(xintercept = median+2*sd,linetype = "dashed", linewidth= 0.5)+
  geom_vline(xintercept = median-2*sd,linetype = "dashed", linewidth = 0.5)+
  scale_x_continuous(breaks=seq(-5,10,5), labels=c(-5,0,5,10))+
  theme_half_open()+
  labs(y="Density",x="Median(Log2PTR)")



####密度散点图绘制####
#准备文件：PRO_n；PTR_n；RNA_n
library(ggplot2)
library(dplyr)
library(viridis) # 使用viridis提供的翠绿色标度：scale_fill_viridis()
library(ggpointdensity) # 绘制密度散点图
library(cowplot) # 图形组合，可以自动对其坐标轴
PRO_n<-PRO[rowSums(is.na(PRO)) <= 4, ]
RNA_n <-RNA[rowSums(is.na(RNA)) <= 4, ]
PRO_n <- apply(PRO_n, c(1, 2), log2)
RNA_n <- apply(RNA_n, c(1, 2), log2)
PRO_n <- as.data.frame(PRO_n)
RNA_n <- as.data.frame(RNA_n)
PRO_n$MAD <- apply(PRO_n[, 1:14], 1, mad, na.rm = TRUE)
RNA_n$MAD <- apply(RNA_n[, 1:14], 1, mad, na.rm = TRUE)
PTR_n$MAD <- apply(PTR_n[, 1:14], 1, mad, na.rm = TRUE)
PRO_n$Median <- apply(PRO_n[, 1:14], 1, median, na.rm = TRUE)
RNA_n$Median <- apply(RNA_n[, 1:14], 1, median, na.rm = TRUE)
PTR_n$Median <- apply(PTR_n[, 1:14], 1, median, na.rm = TRUE)
PRO_P<- PRO_n[,c(15,16)]
RNA_P<- RNA_n[,c(15,16)]
PTR_P<- PTR_n[,c(15,16)]
indices <- c(1734, 3468, 5202, 6936, 8673)
####PRO绘制####
ggplot(data = PRO_P, mapping = aes(x = Median, y = MAD)) +
  geom_pointdensity() +
  scale_color_viridis(option="C",begin = 0.05,end =0.9)+
  guides(colour = guide_colourbar(label = FALSE),colour = guide_colourbar(ticks = FALSE))+
  theme_half_open() + theme(legend.title = element_text(family = "serif",
    colour = NA)) + theme(legend.position = c(0.1, 0.8))+
#  geom_hline(yintercept = PRO_MAD$PRO_MAD[1],color="#CB181D",,linetype = "dashed", linewidth= 0.75)+
#  geom_hline(yintercept = sum(PRO_MAD$PRO_MAD[1:2]),color="#FFFF33",,linetype = "dashed", linewidth= 0.75)+
#  geom_hline(yintercept = sum(PRO_MAD$PRO_MAD[1:3]),color="#238B45",linetype = "dashed", linewidth= 0.75)+
#  geom_hline(yintercept = sum(PRO_MAD$PRO_MAD[1:4]),color="#2171B5",,linetype = "dashed", linewidth= 0.75)+
  theme(legend.key.size = unit(10,"pt"))+
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6))+
  labs(y="MAD protein(log2Expression)",x="Median protein(log2Expression)")
##柱状图绘制
PRO_MAD <-sort(PRO_P$MAD)[indices]
PRO_MAD <-c(PRO_MAD[1], diff(PRO_MAD))
PRO_MAD <-as.data.frame(PRO_MAD)
PRO_MAD$Category <-c("Q1","Q2","Q3","Q4","Q5")
PRO_MAD$Category <-factor(PRO_MAD$Category,levels = rev(c("Q1","Q2","Q3","Q4","Q5"))) 
PRO_MAD$X <-"X"
ggplot(PRO_MAD,aes(y=PRO_MAD,x=X))+
  geom_col(aes(fill=Category),position = "stack",width = 0.2)+
  scale_fill_manual(values = c("#6A51A3","#2171B5","#238B45","#FFFF33","#CB181D"))+
  theme_half_open()+
  theme(legend.position="none")
####RNA绘制####
ggplot(data = RNA_P, mapping = aes(x = Median, y = MAD)) +
  geom_pointdensity(adjust=3) +
  scale_color_viridis(option="C",begin = 0,end =1)+
  guides(colour = guide_colourbar(label = FALSE),colour = guide_colourbar(ticks = FALSE))+
  theme_half_open() + theme(legend.title = element_text(family = "serif",
                                                        colour = NA)) + theme(legend.position = c(0.1, 0.8))+
  theme(legend.key.size = unit(10,"pt"))+
  labs(y="MAD transcript(log2TPM)",x="Median transcript(log2TPM)")
##柱状图绘制
RNA_MAD <-sort(RNA_P$MAD)[indices]
RNA_MAD <-c(RNA_MAD[1], diff(RNA_MAD))
RNA_MAD <-as.data.frame(RNA_MAD)
RNA_MAD$Category <-c("Q1","Q2","Q3","Q4","Q5")
RNA_MAD$Category <-factor(RNA_MAD$Category,levels = rev(c("Q1","Q2","Q3","Q4","Q5"))) 
RNA_MAD$X <-"X"
ggplot(RNA_MAD,aes(y=RNA_MAD,x=X))+
  geom_col(aes(fill=Category),position = "stack",width = 0.2)+
  scale_fill_manual(values = c("#6A51A3","#2171B5","#238B45","#FFFF33","#CB181D"))+
  theme_half_open()+
  theme(legend.position="none")
####PTR绘制####
ggplot(data = PTR_P, mapping = aes(x = Median, y = MAD)) +
  geom_pointdensity(adjust=3) +
  scale_color_viridis(option="C",begin = 0,end =1)+
  guides(colour = guide_colourbar(label = FALSE),colour = guide_colourbar(ticks = FALSE))+
  theme_half_open() + theme(legend.title = element_text(family = "serif",
                                                        colour = NA)) + theme(legend.position = c(0.1, 0.8))+
  theme(legend.key.size = unit(10,"pt"))+
  labs(y="MAD log2PTR",x="Median log2PTR")
##柱状图绘制
PTR_MAD <-sort(PTR_P$MAD)[indices]
PTR_MAD <-c(PTR_MAD[1], diff(PTR_MAD))
PTR_MAD <-as.data.frame(PTR_MAD)
PTR_MAD$Category <-c("Q1","Q2","Q3","Q4","Q5")
PTR_MAD$Category <-factor(PTR_MAD$Category,levels = rev(c("Q1","Q2","Q3","Q4","Q5"))) 
PTR_MAD$X <-"X"
ggplot(PTR_MAD,aes(y=PTR_MAD,x=X))+
  geom_col(aes(fill=Category),position = "stack",width = 0.2)+
  scale_fill_manual(values = c("#6A51A3","#2171B5","#238B45","#FFFF33","#CB181D"))+
  theme_half_open()+
  theme(legend.position="none")




