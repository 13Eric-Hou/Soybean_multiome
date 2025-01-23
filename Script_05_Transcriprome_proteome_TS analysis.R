####TS分数计算####
#加载软件包程序#
rm(list = ls())
library(tidyverse)
library(dplyr)
#加载Tissue-Specificy计算程序
source('D:/Workspace/Glycine_max/Tissue_Specifity/AdaTiSS-main/R/AdaTiSS_fn.R')
#列名准备一下#
col_names <- c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")
#RNA数据导入#
RNA_sample <- read.csv(file = "./Hou/RNA/TPM.csv",row.names = 1)#存在着小于1的数据以及NA值使用0进行代替
temp<-rownames(RNA_sample)
#表达量小于1的单元格视为没有表达量转换为0#
RNA_sample <- data.frame(lapply(RNA_sample, function(x) ifelse(x < 1, 0, x)))
rownames(RNA_sample) <- temp
#PRO数据导入#
PRO_sample <- read.csv(file = "./Hou/PRO/PRO.csv",row.names = 1)
#提取共有基因（12758）#：
RNA_co_quantified_sample <-subset(RNA_sample, rownames(RNA_sample) %in% rownames(PRO_sample))

RNA_co_quantified_sample <- RNA_co_quantified_sample[rowSums(RNA_co_quantified_sample == 0, na.rm = TRUE) != ncol(RNA_co_quantified_sample), ]
PRO_co_quantified_sample <-subset(PRO_sample, rownames(PRO_sample) %in% rownames(RNA_co_quantified_sample))
write.csv(PRO_co_quantified_sample,file = "./Hou/PRO/PRO_co_quantified_sample.csv")
write.csv(RNA_co_quantified_sample,file = "./Hou/RNA/RNA_co_quantified_sample.csv")

rm(RNA_sample)
rm(PRO_sample)
rm(temp)
#PRO数据取对数（matrix）#
PRO_co_quantified_sample <- as.matrix(PRO_co_quantified_sample)
PRO_co_quantified_sample_log = preproc.filter.fn(PRO_co_quantified_sample,dat.type = "intensity", filter.col.prp=1)	#进行对数的计算以及过滤
#RNA数据取对数（matrix）#
RNA_co_quantified_sample <- as.matrix(RNA_co_quantified_sample)
RNA_co_quantified_sample_log = preproc.filter.fn(RNA_co_quantified_sample,dat.type = "TPM or RPKM", proc.zero ="perturbed by a small value",exp.thres=0)
#准备组织水平的文件#
Sample_info <- read.csv("./Hou/Sample_info_Sample.csv")
PRO_tiss.abd <- tiss.abd.fn(PRO_co_quantified_sample_log, Sample_info)
RNA_tiss.abd <- tiss.abd.fn(RNA_co_quantified_sample_log, Sample_info)
PRO_tiss.abd <- PRO_tiss.abd[, c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")]
RNA_tiss.abd <- RNA_tiss.abd[, c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")]
#此处是三个数值取中位数，两个数值取平均数，一个数值就是这个数值#（尝试三个数值也是求平均值）
for (i in 1:14) 
  {
    for (j in 1:nrow(PRO_co_quantified_sample_log)) {
      Temp <- c(PRO_co_quantified_sample_log[j,3*i-1],PRO_co_quantified_sample_log[j,3*i-2],PRO_co_quantified_sample_log[j,3*i])
      if (sum(is.na(Temp))==3) 
      {
        PRO_tiss.abd[j,i] <- NA
      }
      if (sum(is.na(Temp))!=3) 
      {
        PRO_tiss.abd[j,i] <- mean(Temp,na.rm = T)
      }
    }
  }
for (i in 1:14) 
{
  for (j in 1:nrow(RNA_co_quantified_sample_log)) {
    Temp <- c(RNA_co_quantified_sample_log[j,3*i-1],RNA_co_quantified_sample_log[j,3*i-2],RNA_co_quantified_sample_log[j,3*i])
    if (sum(is.na(Temp))==3) 
    {
      RNA_tiss.abd[j,i] <- NA
    }
    if (sum(is.na(Temp))!=3) 
    {
       RNA_tiss.abd[j,i] <- mean(Temp,na.rm = T)
    }
  }
}





#计算TS分数#
PRO_result = AdaTiSS(PRO_co_quantified_sample_log,tiss.abd = PRO_tiss.abd)
RNA_result = AdaTiSS(RNA_co_quantified_sample_log,tiss.abd = RNA_tiss.abd,adjust.opt = 0)
#组织TS分数提取#
PRO_Score_Sample <- PRO_result$ada.s
RNA_Score_Sample <- RNA_result$ada.s
PRO_Score_Tissue <- PRO_result$ada.z
RNA_Score_Tissue_raw <- RNA_result$ada.z
RNA_Score_Tissue <- RNA_result$ada.z
#此处是三个数值取中位数，两个数值取平均数，一个数值就是这个数值#（尝试三个数值也是求平均值）#
for (i in 1:14) 
{
  for (j in 1:nrow(PRO_Score_Sample)) {
    Temp <- c(PRO_Score_Sample[j,3*i-1],PRO_Score_Sample[j,3*i-2],PRO_Score_Sample[j,3*i])
    if (sum(is.na(Temp))==3) 
    {
      PRO_Score_Tissue[j,i] <- NA
    }
    if (sum(is.na(Temp))!=3) 
    {
      PRO_Score_Tissue[j,i] <- mean(Temp,na.rm = T)
    }
  }
}
for (i in 1:14) 
{
  for (j in 1:nrow(RNA_Score_Sample)) {
    Temp <- c(RNA_Score_Sample[j,3*i-1],RNA_Score_Sample[j,3*i-2],RNA_Score_Sample[j,3*i])
    if (sum(is.na(Temp))==3) 
    {
      RNA_Score_Tissue_raw[j,i] <- NA
    }
    if (sum(is.na(Temp))!=3) 
    {
      RNA_Score_Tissue_raw[j,i] <- mean(Temp,na.rm = T)
    }
  }
}
for (i in 1:14) 
{
  for (j in 1:nrow(RNA_Score_Sample)) {
    Temp <- c(RNA_Score_Sample[j,3*i-1],RNA_Score_Sample[j,3*i-2],RNA_Score_Sample[j,3*i])
    if (sum(is.na(Temp))==3) 
    {
      RNA_Score_Tissue[j,i] <- NA
    }
    if (sum(is.na(Temp))!=3) 
    {
      RNA_Score_Tissue[j,i] <- mean(Temp,na.rm = T)
    }
  }
}
#把这三个表格大于10的单元格转换为10#
PRO_Score_Tissue <- apply(PRO_Score_Tissue, 2, function(x) ifelse(x > 10, 10, x))
RNA_Score_Tissue <- apply(RNA_Score_Tissue, 2, function(x) ifelse(x > 10, 10, x))
RNA_Score_Tissue_raw <- apply(RNA_Score_Tissue_raw, 2, function(x) ifelse(x > 10, 10, x))

write.csv(PRO_Score_Sample,file = "./Hou/PRO/PRO_Score_Sample.csv")
write.csv(RNA_Score_Sample,file = "./Hou/RNA/RNA_Score_Sample.csv")

PRO_Score_Tissue <- PRO_Score_Tissue[, c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")]
RNA_Score_Tissue <- RNA_Score_Tissue[, c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")]
rm(PRO_tiss.abd)
rm(RNA_tiss.abd)
rm(Sample_info)
rm(PRO_co_quantified_sample_log)
rm(RNA_co_quantified_sample_log)
write.csv(PRO_Score_Tissue,file = "./Hou/PRO/PRO_Score_tissue.csv")
####TS分数过滤####
##RNA过滤1:组织水平的TPM值都小于1，设置为all...、
#平均到组织文件（平均值）
RNA_raw_tiss_abd <- data.frame(matrix(nrow = 12758, ncol = 14))
for (i in 1:14) 
{
  for (j in 1:nrow(RNA_co_quantified_sample)) {
    Temp <- c(RNA_co_quantified_sample[j,3*i-1],RNA_co_quantified_sample[j,3*i-2],RNA_co_quantified_sample[j,3*i])
    if (sum(is.na(Temp))==3) 
    {
      RNA_raw_tiss_abd [j,i] <- NA
    }
    if (sum(is.na(Temp))!=3) 
    {
      RNA_raw_tiss_abd [j,i] <- mean(Temp,na.rm = T)
    }
  }
}
RNA_raw_tiss_abd <- apply(RNA_raw_tiss_abd, 2, function(x) ifelse(x == 0, NA, x))
PRO_raw_tiss_abd <- data.frame(matrix(nrow = 12758, ncol = 14))
for (i in 1:14) 
{
  for (j in 1:nrow(PRO_co_quantified_sample)) {
    Temp <- c(PRO_co_quantified_sample[j,3*i-1],PRO_co_quantified_sample[j,3*i-2],PRO_co_quantified_sample[j,3*i])
    if (sum(is.na(Temp))==3) 
    {
      PRO_raw_tiss_abd [j,i] <- NA
    }
    if (sum(is.na(Temp))!=3) 
    {
      PRO_raw_tiss_abd [j,i] <- mean(Temp,na.rm = T)
    }
  }
}
PRO_raw_tiss_abd <- apply(PRO_raw_tiss_abd, 2, function(x) ifelse(x == 0, NA, x))
#给组织水平的RNA文件行名以及列名
row.names(RNA_raw_tiss_abd) <- row.names(RNA_co_quantified_sample)
colnames(RNA_raw_tiss_abd) <- col_names
row.names(PRO_raw_tiss_abd) <- row.names(PRO_co_quantified_sample)
colnames(PRO_raw_tiss_abd) <- col_names
#在所有组织的平均表达量小于1的筛选出来,并把RNA-score——tissue文件替换
RNA_raw_tiss_abd <- apply(RNA_raw_tiss_abd, 2, function(x) ifelse(is.na(x), 0, x))
RNA_filter_1 <- RNA_raw_tiss_abd[rowSums(RNA_raw_tiss_abd< 1) == ncol(RNA_raw_tiss_abd), ]
LIST <- rownames(RNA_filter_1)
RNA_Score_Tissue[LIST,] <- "NA_all_tissues_tpm_less_1"
rm(RNA_filter_1)
##RNA过滤2：结合阈值进行过滤
#表达量判断：表达量小于1的该单元格设为该列的列名，大于1则视为NA
Expression_adjust <- data.frame(matrix("", nrow = nrow(RNA_raw_tiss_abd), ncol = ncol(RNA_raw_tiss_abd)))
rownames(Expression_adjust) <- rownames(RNA_raw_tiss_abd)
colnames(Expression_adjust) <- colnames(RNA_raw_tiss_abd)
for (i in 1:nrow(RNA_raw_tiss_abd)) {
  for (j in 1:ncol(RNA_raw_tiss_abd)) {
    if (RNA_raw_tiss_abd[i, j] >= 1) {
      Expression_adjust[i, j] <- NA
    } else {
      Expression_adjust[i, j] <- colnames(RNA_raw_tiss_abd)[j]
    }
  }
}

#TS分数过滤：如果某个单元格的TS分数大于2.5则改为这一列的列名，否则转换为NA
TS_adjust <- data.frame(matrix("", nrow = nrow(RNA_Score_Tissue), ncol = ncol(RNA_Score_Tissue)))
rownames(TS_adjust) <- rownames(RNA_Score_Tissue)
colnames(TS_adjust) <- colnames(RNA_Score_Tissue)
RNA_Score_Tissue_NA <- RNA_Score_Tissue 
RNA_Score_Tissue_NA[RNA_Score_Tissue_NA == "NA_all_tissues_tpm_less_1"] <- -666
for (i in 1:nrow(RNA_Score_Tissue_NA)) {
  for (j in 1:ncol(RNA_Score_Tissue_NA)) {
    if (RNA_Score_Tissue_NA[i, j] >=2.5) {
      TS_adjust[i, j] <- colnames(RNA_Score_Tissue)[j]
    } else {
      TS_adjust[i, j] <- NA
    }
  }
}
#双判断:判断TS分数大于2.5但是原始表达量小于1
RNA_filtered_2 <- data.frame(matrix(NA, nrow = nrow(RNA_Score_Tissue), ncol = ncol(RNA_Score_Tissue)))
rownames(RNA_filtered_2) <- rownames(RNA_Score_Tissue)
colnames(RNA_filtered_2) <- colnames(RNA_Score_Tissue)
for (i in 1:nrow(RNA_filtered_2)) {
  for (j in 1:ncol(RNA_filtered_2)) {
    if (!is.na(TS_adjust[i, j]) && !is.na(Expression_adjust[i, j]) && TS_adjust[i, j] == Expression_adjust[i, j]) {
      RNA_filtered_2[i, j] <- ";NA_raw_tpm_less_1"
    }
  }
}
#追加
RNA_filtered_2[is.na(RNA_filtered_2)] <-"" #NA转换为空，非空单元格，空单元格追加会追加上NA
for (i in 1:nrow(RNA_filtered_2)) 
{
  for (j in 1:ncol(RNA_filtered_2)) 
  {
    RNA_Score_Tissue [i, j] <- paste(RNA_Score_Tissue[i, j], RNA_filtered_2[i, j], sep = "")
  }
}
rm(RNA_filtered_2)
rm(Expression_adjust)
rm(TS_adjust)
rm(RNA_Score_Tissue_NA)
#write.csv(RNA_Score_Tissue,file = "./Hou/RNA/RNA_Score_Tissue.csv")
####分类####
####蛋白质水平分类####
PRO_filtered <- PRO_Score_Tissue
#PRO_filtered[grepl("NA", PRO_filtered)] <- NA
PRO_filtered <- as.data.frame(PRO_filtered)
##PRO候选HousekeepinG蛋白的筛选：至少要在一个组织中鉴定到一次（补缺前）
PRO_housekeeping_raw <- PRO_raw_tiss_abd[rowSums(is.na(PRO_raw_tiss_abd)) == 0, ]
housekeeping_raw_list<- row.names(PRO_housekeeping_raw)
PRO_housekeeping_raw <-subset(PRO_filtered, rownames(PRO_filtered) %in% housekeeping_raw_list) 
PRO_housekeeping_raw <- as.data.frame(PRO_housekeeping_raw)
PRO_housekeeping <-PRO_housekeeping_raw[rowSums(is.na(PRO_housekeeping_raw)) == 0 & apply(PRO_housekeeping_raw, 1, function(row) all(row < 2)), ]
#删除中间过渡文件
rm(PRO_housekeeping_raw)
#write.csv(PRO_housekeeping,file = "./Hou/PRO/PRO_housekeeping.csv")
##Enriched & Specific
PRO_filtered <- as.data.frame(PRO_filtered)
PRO_filtered$Specific <- apply(PRO_filtered, 1, function(row) sum(row >=4, na.rm = TRUE))
PRO_filtered$Enriched <- apply(PRO_filtered, 1, function(row) sum(row >= 2.5, na.rm = TRUE))
PRO_Specific <- PRO_filtered[PRO_filtered[, 15] > 0& PRO_filtered[, 15]==PRO_filtered[,16], ]
PRO_Enriched <- PRO_filtered[PRO_filtered[, 16] >0,]
PRO_Enriched_not_specific <- PRO_Enriched[-which(rownames(PRO_Enriched) %in% rownames(PRO_Specific )), ]
PRO_Specific<-PRO_Specific[,-15]
PRO_Specific<-PRO_Specific[,-15]
PRO_Enriched<-PRO_Enriched[,-15]
PRO_Enriched<-PRO_Enriched[,-15]
PRO_Enriched_not_specific<-PRO_Enriched_not_specific[,-15]
PRO_Enriched_not_specific<-PRO_Enriched_not_specific[,-15]
#write.csv(PRO_Specific,file = "./Hou/PRO/PRO_Specific.csv")
#write.csv(PRO_Enriched_not_specific,file = "./Hou/PRO/PRO_Enriched_not_specific.csv")
##others
PRO_others <- PRO_filtered[-which(rownames(PRO_filtered) %in% rownames(PRO_housekeeping)), ]
PRO_others <- PRO_others[-which(rownames(PRO_others) %in% rownames(PRO_Enriched)), ]
PRO_others<-PRO_others[,-15]
PRO_others<-PRO_others[,-15]
#write.csv(PRO_others,file = "./Hou/PRO/PRO_others.csv")
rm(PRO_filtered)
####RNA水平分类####
RNA_filtered <- RNA_Score_Tissue
#把所有含字符NA的单元格转换为NA
RNA_filtered[grepl("NA", RNA_filtered)] <- NA
RNA_filtered <- as.data.frame(RNA_filtered)
for (col in 1:ncol(RNA_filtered)) {
  # 将每一列转换为数值型
  RNA_filtered[, col] <- as.numeric(RNA_filtered[, col])
}
#RNA Housekeeping判断：
RNA_housekeeping_raw <- RNA_raw_tiss_abd[rowSums(RNA_raw_tiss_abd == 0) == 0, ]
housekeeping_raw_list<- row.names(RNA_housekeeping_raw)
RNA_housekeeping_raw <-subset(RNA_filtered, rownames(RNA_filtered) %in% housekeeping_raw_list) 
RNA_housekeeping_raw <- as.data.frame(RNA_housekeeping_raw)
RNA_housekeeping <-RNA_housekeeping_raw[rowSums(is.na(RNA_housekeeping_raw)) == 0 & apply(RNA_housekeeping_raw, 1, function(row) all(row < 2)), ]
##Enriched & Specific
RNA_filtered <- as.data.frame(RNA_filtered)
RNA_filtered$Specific <- apply(RNA_filtered, 1, function(row) sum(row >4, na.rm = TRUE))
RNA_filtered$Enriched<- apply(RNA_filtered, 1, function(row) sum(row > 2.5, na.rm = TRUE))
RNA_Specific <- RNA_filtered[RNA_filtered[, 15] > 0&RNA_filtered[, 15]==RNA_filtered[, 16], ]
RNA_Enriched <- RNA_filtered[RNA_filtered[, 16] >0,]
RNA_Enriched_not_specific <- RNA_Enriched[-which(rownames(RNA_Enriched) %in% rownames(RNA_Specific )), ]
RNA_Specific<-RNA_Specific[,-15]
RNA_Specific<-RNA_Specific[,-15]
RNA_Enriched<-RNA_Enriched[,-15]
RNA_Enriched<-RNA_Enriched[,-15]
RNA_Enriched_not_specific<-RNA_Enriched_not_specific[,-15]
RNA_Enriched_not_specific<-RNA_Enriched_not_specific[,-15]
##others#
RNA_others <- RNA_filtered[-which(rownames(RNA_filtered) %in% rownames(RNA_housekeeping )), ]
RNA_others <-RNA_others[-which(rownames(RNA_others) %in% rownames(RNA_Enriched)), ]
RNA_others<-RNA_others[,-15]
RNA_others<-RNA_others[,-15]
rm(RNA_housekeeping_raw)
rm(RNA_filtered)
rm(housekeeping_raw_list)
rm(RNA_raw_tiss_abd,RNA_co_quantified_sample)
rm(RNA_Score_Tissue_raw)
#导出文件#
write.csv(RNA_Enriched,file = "./Hou/RNA/RNA_Enriched.csv")
write.csv(RNA_housekeeping,file = "./Hou/RNA/RNA_housekeeping.csv")
write.csv(RNA_Enriched_not_specific,file = "./Hou/RNA/RNA_Enriched_not_specific.CSV")
write.csv(RNA_Specific,file = "D:/Workspace/Glycine_max/Tissue_Specifity/Hou/RNA/RNA_Specific.csv")
write.csv(RNA_others,file = "./Hou/RNA/RNA_others.csv")