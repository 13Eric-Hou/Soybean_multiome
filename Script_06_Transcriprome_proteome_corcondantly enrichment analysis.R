##PRO与RNA的一致与不一致的问题
##准备文件：PRO_Score_Tissue与RNA_Score_Tissue
#both enriched:某个蛋白他的TS得分都大于1.25
#RNA enriched:RNA的蛋白得分大于1.25，而且RNA的TS得分比PRO得分高0.75
#PRO enriched:PRO的蛋白得分大于1.25，而且PRO的TS得分比RNA得分高0.75
library(reshape2)
library(ggplot2)
library(RColorBrewer)
rm(list=ls())
#导入文件#
RNA_Score_Tissue <- read.csv(file = "../Tissue_Specifity/Hou/RNA/RNA_Score_Tissue.csv",row.names = 1)
PRO_Score_Tissue <- read.csv(file = "../Tissue_Specifity/Hou/PRO/PRO_Score_Tissue.csv",row.names = 1)
Concordance_adjust <- data.frame(matrix("", nrow = nrow(RNA_Score_Tissue), ncol = ncol(RNA_Score_Tissue)))
#把含有NA的单元格转换为NA
RNA_CON <- RNA_Score_Tissue
RNA_CON <- as.matrix(RNA_CON)
RNA_CON[grepl("NA", RNA_CON)] <- NA
RNA_CON <- as.data.frame(RNA_CON)
#把含有NA的单元格转换为NA
PRO_CON <- PRO_Score_Tissue

for (col in 1:ncol(RNA_CON)) {
  # 将每一列转换为数值型
  RNA_CON[, col] <- as.numeric(RNA_CON[, col])
}
for (col in 1:ncol(PRO_CON)) {
  # 将每一列转换为数值型
  PRO_CON[, col] <- as.numeric(PRO_CON[, col])
}

#判断文件命名
rownames(Concordance_adjust) <- rownames(RNA_Score_Tissue)
colnames(Concordance_adjust) <- colnames(RNA_Score_Tissue)

for (i in 1:nrow(PRO_CON)) {
  for (j in 1:ncol(PRO_CON)) {
    # 检查两个文件相同位置的单元格是否都大于1.25
    if (PRO_CON[i, j] > 2.5 && RNA_CON[i, j] > 2.5 &&!is.na(PRO_CON[i, j])&&!is.na(RNA_CON[i, j])) {
      # 赋值为"both_enriched"
      Concordance_adjust[i, j] <- "both_enriched"
    }
  }
}


for (i in 1:nrow(PRO_CON)) {
  for (j in 1:ncol(PRO_CON)) {
    # 检查文件1的单元格是否大于1.25，且比文件2的单元格大0.75
    if (!is.na(PRO_CON[i, j]) &&PRO_CON[i, j] > 2.5&& RNA_CON[i, j] < 2.5&&!is.na(RNA_CON[i, j]) && (PRO_CON[i, j] >= RNA_CON[i, j] + 1.5)) {
      # 赋值为"PRO_enriched"
      Concordance_adjust[i, j] <- "PRO_enriched"
    }
    
    # 检查文件2的单元格是否大于1.25，且比文件1的单元格大0.75
    if (!is.na(RNA_CON[i, j]) && RNA_CON[i, j] > 2.5&& PRO_CON[i, j]< 2.5 && !is.na(PRO_CON[i, j]) && (RNA_CON[i, j] >= PRO_CON[i, j] + 1.5)) {
      # 赋值为"RNA_enriched"
      Concordance_adjust[i, j] <- "RNA_enriched"
    }
  }
}
Concordance_result <- data.frame(matrix("", nrow = 14, ncol =3))
rownames(Concordance_result) <- colnames(RNA_Score_Tissue)
colnames(Concordance_result) <- c("PRO_enriched","RNA_enriched","Both_enriched")
for (i in 1:14) {
  # 统计"PRO_enriched"单元格数量
  Concordance_result[i, 1] <- sum(Concordance_adjust[, i] == "PRO_enriched", na.rm = TRUE)
  
  # 统计"RNA_enriched"单元格数量
  Concordance_result[i, 2] <- sum(Concordance_adjust[, i] == "RNA_enriched", na.rm = TRUE)
  
  # 统计"both_enriched"单元格数量
  Concordance_result[i, 3] <- sum(Concordance_adjust[, i] == "both_enriched", na.rm = TRUE)
}
Concordance_result$Tissue <- rownames(Concordance_result)
Concordance_result<- melt(Concordance_result,id.vars='Tissue')
Concordance_result$value <- as.numeric(Concordance_result$value)
Concordance_result$Tissue <- factor(Concordance_result$Tissue,levels =c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
                                                                        "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2"))
##绘图
ggplot(data=Concordance_result,aes(Tissue,value,fill=variable))+
  geom_col(position="stack", color="black", width=0.8,linewidth=0)+
  scale_fill_manual(values = c("#74C476","#6BAED6","#9E9AC8"))+
  theme(axis.text = element_text(hjust = 1))+
  theme_test() + 
  theme(axis.text.x = element_text(angle = 45)) + theme(legend.position = c(0.8, 0.95))+labs(y = "Count") +
  theme(axis.text.x = element_text(vjust = 0.7))+ 
  labs(y="Number", x="Tissue") +
  theme(legend.position = c(0.68, 0.96),legend.key.size = unit(7, "pt"))+ 
  theme(legend.text = element_text(size = 7),
        legend.title = element_text(colour = NA))+ 
  theme(legend.background = element_rect(fill = NA))+ 
  theme(axis.text.x = element_text(size = 7))
