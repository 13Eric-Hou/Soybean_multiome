####每个组织表达量Rank排列，搭配Y轴累计表达量####
rm(list = ls())
options(scipen = 100000)
#导入表达量数据整理，添加evidence列#
PRO <- read.csv("D:/Workspace/Glycine_max/Data/PRO_Raw_Gene_Tissue.csv",row.names = 1)
RNA <- read.csv("D:/Workspace/Glycine_max/Data/RNA_Tissue.csv",row.names = 1)
PRO$omics_level <-"PRO"
PRO$Gene_name <- rownames(PRO)
rownames(PRO) <- 1:nrow(PRO)
RNA$omics_level <-"RNA"
RNA$Gene_name <- rownames(RNA)
rownames(RNA) <- 1:nrow(RNA)

####准备每个组织的文件####
Data_1_PRO <- PRO[!is.na(PRO[, 1]), c(1, 15, 16)]
Data_2_PRO <- PRO[!is.na(PRO[, 2]), c(2, 15, 16)]
Data_3_PRO <- PRO[!is.na(PRO[, 3]), c(3, 15, 16)]
Data_4_PRO <- PRO[!is.na(PRO[, 4]), c(4, 15, 16)]
Data_5_PRO <- PRO[!is.na(PRO[, 5]), c(5, 15, 16)]
Data_6_PRO <- PRO[!is.na(PRO[, 6]), c(6, 15, 16)]
Data_7_PRO <- PRO[!is.na(PRO[, 7]), c(7, 15, 16)]
Data_8_PRO <- PRO[!is.na(PRO[, 8]), c(8, 15, 16)]
Data_9_PRO <- PRO[!is.na(PRO[, 9]), c(9, 15, 16)]
Data_10_PRO <-PRO[!is.na(PRO[, 10]), c(10, 15, 16)]
Data_11_PRO <- PRO[!is.na(PRO[, 11]), c(11, 15, 16)]
Data_12_PRO <- PRO[!is.na(PRO[, 12]), c(12, 15, 16)]
Data_13_PRO <- PRO[!is.na(PRO[, 13]), c(13, 15, 16)]
Data_14_PRO <- PRO[!is.na(PRO[, 14]), c(14, 15, 16)]
Data_1_RNA <- RNA[!is.na(RNA[, 1]), c(1, 15, 16)]
Data_2_RNA <- RNA[!is.na(RNA[, 2]), c(2, 15, 16)]
Data_3_RNA <- RNA[!is.na(RNA[, 3]), c(3, 15, 16)]
Data_4_RNA <- RNA[!is.na(RNA[, 4]), c(4, 15, 16)]
Data_5_RNA <- RNA[!is.na(RNA[, 5]), c(5, 15, 16)]
Data_6_RNA <- RNA[!is.na(RNA[, 6]), c(6, 15, 16)]
Data_7_RNA <- RNA[!is.na(RNA[, 7]), c(7, 15, 16)]
Data_8_RNA <- RNA[!is.na(RNA[, 8]), c(8, 15, 16)]
Data_9_RNA <- RNA[!is.na(RNA[, 9]), c(9, 15, 16)]
Data_10_RNA <-RNA[!is.na(RNA[, 10]), c(10, 15, 16)]
Data_11_RNA <- RNA[!is.na(RNA[, 11]), c(11, 15, 16)]
Data_12_RNA <- RNA[!is.na(RNA[, 12]), c(12, 15, 16)]
Data_13_RNA <- RNA[!is.na(RNA[, 13]), c(13, 15, 16)]
Data_14_RNA <- RNA[!is.na(RNA[, 14]), c(14, 15, 16)]
#降序
for (i in 1:14) {
  assign(paste0("Data_", i, "_PRO"), get(paste0("Data_", i, "_PRO"))[order(-get(paste0("Data_", i, "_PRO"))[, 1]), ])
}
for (i in 1:14) {
  assign(paste0("Data_", i, "_RNA"), get(paste0("Data_", i, "_RNA"))[order(-get(paste0("Data_", i, "_RNA"))[, 1]), ])
}
#加Rank
# 为每个数据框添加一列Rank，内容为1到文件行数的连续数字
# 为每个数据框添加一列Rank，内容为1到文件行数的连续数字
for (i in 1:14) {
  rows <- nrow(get(paste0("Data_", i, "_PRO")))
  assign(paste0("Data_", i, "_PRO"), cbind(get(paste0("Data_", i, "_PRO")), Rank = seq(1, rows)))
}
for (i in 1:14) {
  rows <- nrow(get(paste0("Data_", i, "_RNA")))
  assign(paste0("Data_", i, "_RNA"), cbind(get(paste0("Data_", i, "_RNA")), Rank = seq(1, rows)))
}
# 为每个数据框添加一列Intensity[%]，内容为NA
for (i in 1:14) {
  data <- get(paste0("Data_", i, "_PRO"))
  data$Intensity <- (data[, 1] / sum(data[, 1])) * 100
  assign(paste0("Data_", i, "_PRO"), data)
}
for (i in 1:14) {
  data <- get(paste0("Data_", i, "_RNA"))
  data$Intensity <- (data[, 1] / sum(data[, 1])) * 100
  assign(paste0("Data_", i, "_RNA"), data)
}
for (i in 1:14) {
  data <- get(paste0("Data_", i, "_PRO"))
  data$SUM <- cumsum(data[, 5])
  assign(paste0("Data_", i, "_PRO"), data)
}
for (i in 1:14) {
  data <- get(paste0("Data_", i, "_RNA"))
  data$SUM <- cumsum(data[, 5])
  assign(paste0("Data_", i, "_RNA"), data)
}
# 循环遍历每个数据框，将PRO和RNA数据框合并存放在新的数据框中
for (i in 1:14) {
  data_PRO <- get(paste0("Data_", i, "_PRO"))
  data_RNA <- get(paste0("Data_", i, "_RNA"))
  assign(paste0("Data_", i), rbind(data_PRO, data_RNA))
}
# 只保留所有蛋白基因层次#RNA基因有删除
# 循环遍历每个文件
for (i in 1:14) {
  # 获取文件名
  file_name <- paste0("Data_", i)
  # 读取文件
  data <- get(file_name)
  # 计算第四列的最大值
  max_value <- max(data[data[, 2] == "PRO", ][,4])
  # 删除满足条件且第四列等于最大值的行
  condition <- data[, 2] == "RNA" & data[, 4] > max_value
  data <- data[!condition,]
  # 更新文件
  assign(file_name, data)
}
####绘图####
library(ggplot2)
library(tidyverse)
library(ggthemes)
library(RColorBrewer)
for (i in 1:14)
  {
  # 获取文件名
  file_name <- paste0("Data_", i)
  # 读取文件
  data <- get(file_name)
  data$Rank <- log10(data$Rank)
  # 更新文件
  assign(file_name, data)
}
###绘图###
for (i in 1:14) {
  file_name <- paste0("Data_", i)
  data <- get(file_name)
  if (i==11|2|3|6) {
    temp <- ggplot(data,aes(x=Rank,y=SUM))+
      geom_point(aes(color=omics_level),size=0.6)+
      scale_color_manual(values = c("#74C476","#6BAED6"))+
      scale_y_continuous(breaks=seq(0,100,20), labels=c(0,20,40,60,80,100))+
      scale_x_continuous(breaks=seq(0,4,1), labels=c("1E0","1E1","1E2","1E3","1E4"))+
      labs(y="Cumulative abundance[%]")+
      theme_test() +
      
      theme(legend.title = element_text(size = 0)) + theme(legend.background = element_rect(fill = NA), legend.key.width = unit(15, "pt")) + 
      theme(legend.text = element_text(size = 10)) + theme(legend.position = "none")
    file_name <- paste0("P_", i)
    
    assign(file_name, temp)
  }
  if (i!=11|2|3|6) {
    temp <- ggplot(data,aes(x=Rank,y=SUM))+
      geom_point(aes(color=omics_level),size=0.6)+
      scale_color_manual(values = c("#74C476","#6BAED6"))+
      scale_y_continuous(breaks=seq(0,100,20), labels=c(0,20,40,60,80,100))+
      scale_x_continuous(breaks=seq(0,4,1), labels=c("1E0","1E1","1E2","1E3","1E4"))+
      labs(y=NULL)+
      theme_test() +
      
      theme(legend.title = element_text(size = 0)) + theme(legend.background = element_rect(fill = NA), legend.key.width = unit(15, "pt")) + 
      theme(legend.text = element_text(size = 10)) + theme(legend.position = "none")
    file_name <- paste0("P_", i)
    
    assign(file_name, temp)
  }
 
}
library(patchwork)
(P_11|P_13|P_14|P_1)/(P_2|P_7|P_4|P_8)/(P_3|P_6|P_10|P_12)
P_5|P_9




####基因整理####
Tidy <- data.frame()
for (i in 1:14) {
  # 获取文件名
  file_name <- paste0("Data_", i)
  # 读取文件
  data <- get(file_name)
  top_10_rows <- head(data[order(data[, 4]), ], 10)
  top_10_rows <- top_10_rows[,-1]
  top_10_rows$Tissue <- i
  Tidy<- rbind(Tidy, top_10_rows)
}
write.csv(Tidy,file = "./Tidy.csv")



####TOP100####
#把Rank变化为1，2，3，4等###
for (i in 1:14) {
  # 获取文件名
  file_name <- paste0("Data_", i)
  # 读取文件
  data <- get(file_name)
  data$Rank <- 10^data$Rank
  assign(paste0("Data_", i), data)
}
Tissue <- c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")
Category <- c("Shared","Not Shared")
Top100 <-data.frame(Tissue = rep(Tissue, times = 2),Category = rep(Category, each = 14),count=NA)
for (i in 1:14) 
  {
  file_name <- paste0("Data_", i)
  data <- get(file_name)
  data <-  data[data$Rank >= 1 & data$Rank <= 100, ]
  Top100[i,3] <- (200-length(unique(data$Gene_name)))
  Top100[i+14,3] <- (length(unique(data$Gene_name))-100)
}
Top100$count <- as.numeric(Top100$count)
Top100$Category <- factor(Top100$Category,level=c("Shared","Not Shared"))
Top100$Tissue <- factor(Top100$Tissue ,levels =rev(c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                                 "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2")))
##画图##
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggthemes)
library(ggcor)
library(cowplot)
ggplot(data=Top100,aes(Tissue,count,fill=Category))+
  geom_col(position="stack", color=NA, width=0.9,linewidth=0)+
  theme_half_open()+
  scale_fill_manual(values = c("#41AB5D","#FDB462"))+
  coord_flip()+ 
  theme(axis.text.y = element_text(size = 9))+
  labs(x="Tissue") + 
  labs(y="TOP 100 Count")+
  theme( legend.key.height = unit(10, "pt"),legend.key.width = unit(10, "pt") ) + 
  theme(legend.text = element_text(size = 10),legend.title = element_text(size = 10))+
  theme(legend.position = "bottom", legend.direction = "horizontal")


####图3：不同组织前丰度基因累计柱状图排列####
#寻找到每个组织中表达量最高的蛋白质#

TOP <- data.frame()
for (i in 1:14) 
{
  file_name <- paste0("Data_", i)
  data <- get(file_name)
  colnames(data)[1] <- "Expression"
  temp <- data[data$omics_level == "PRO" & data$Rank == 1, ]
  TOP <- rbind(TOP,temp)
}
Gene <-TOP$Gene_name
Gene <- unique(Gene)
TOP[TOP$Gene_name== "SoyC11_04G086000", "Gene_name"] <- "PSII D2"
TOP[TOP$Gene_name== "SoyC11_16G065900", "Gene_name"] <- "P91"
TOP[TOP$Gene_name== "SoyC11_17G028100", "Gene_name"] <- "SAM22"
TOP[TOP$Gene_name== "SoyC11_10G187900", "Gene_name"] <- "FRK5"
TOP[TOP$Gene_name== "SoyC11_10G207800", "Gene_name"] <- "LB"
TOP[TOP$Gene_name== "SoyC11_07G065100", "Gene_name"] <- "MEF31"
TOP[TOP$Gene_name== "SoyC11_06G058100", "Gene_name"] <- "CesA8"
TOP7 <- data.frame()
for (i in 1:14) 
{
  file_name <- paste0("Data_", i)
  data <- get(file_name)
  colnames(data)[1] <- "Expression"
  data$Tissue <- Tissue[i]
  temp <- data[data$Gene_name %in% Gene&data[,2]=="PRO", ]
  TOP7 <- rbind(TOP7,temp)
}
TOP7 <- TOP7[,-c(1,2,4,6)]
####
TOP7[TOP7$Gene_name== "SoyC11_04G086000", "Gene_name"] <- "PSII D2"
TOP7[TOP7$Gene_name== "SoyC11_16G065900", "Gene_name"] <- "P91"
TOP7[TOP7$Gene_name== "SoyC11_17G028100", "Gene_name"] <- "SAM22"
TOP7[TOP7$Gene_name== "SoyC11_10G187900", "Gene_name"] <- "FRK5"
TOP7[TOP7$Gene_name== "SoyC11_10G207800", "Gene_name"] <- "LB"
TOP7[TOP7$Gene_name== "SoyC11_07G065100", "Gene_name"] <- "MEF31"
TOP7[TOP7$Gene_name== "SoyC11_06G058100", "Gene_name"] <- "CesA8"
library(ggsci)
TOP7$Gene_name <- factor(TOP7$Gene_name,levels = c("PSII D2","P91","FRK5","LB","MEF31","CesA8","SAM22"))

TOP7$Tissue <- factor(TOP7$Tissue,levels = rev(c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                                 "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2")))
ggplot(TOP7,aes(Tissue,Intensity,fill=Gene_name))+
  geom_col(position="stack", color=NA, width=0.9,linewidth=0)+
  theme_half_open()+
  scale_fill_manual(values = c("#41AB5D","#FDB462","#807DBA" ,"#B15928","#EF3B2C","#DD3497","#4292C6"))+
  coord_flip()+ 
  theme(axis.text.y = element_text(size = 9))+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.key.size = unit(8, "pt"))+
  labs(x="Tissue")+ 
  guides(fill = guide_legend(title = "Top 1 Proteins"))+
  theme(legend.text = element_text(size = 8))+ theme(legend.title = element_text(size = 8)) + theme(legend.position = c(0.72, 0.9))
