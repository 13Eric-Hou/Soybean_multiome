####TS SCORE ANALYSIS####

rm(list = ls())
library(tidyverse)
library(dplyr)

source('D:/Workspace/Glycine_max/Tissue_Specifity/AdaTiSS-main/R/AdaTiSS_fn.R')

col_names <- c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")
#LOAD data#
RNA_sample <- read.csv(file = "./Hou/RNA/TPM.csv",row.names = 1)#存在着小于1的数据以及NA值使用0进行代替
temp<-rownames(RNA_sample)

RNA_sample <- data.frame(lapply(RNA_sample, function(x) ifelse(x < 1, 0, x)))
rownames(RNA_sample) <- temp

PRO_sample <- read.csv(file = "./Hou/PRO/PRO.csv",row.names = 1)
#co_quantified genes#：
RNA_co_quantified_sample <-subset(RNA_sample, rownames(RNA_sample) %in% rownames(PRO_sample))

RNA_co_quantified_sample <- RNA_co_quantified_sample[rowSums(RNA_co_quantified_sample == 0, na.rm = TRUE) != ncol(RNA_co_quantified_sample), ]
PRO_co_quantified_sample <-subset(PRO_sample, rownames(PRO_sample) %in% rownames(RNA_co_quantified_sample))
write.csv(PRO_co_quantified_sample,file = "./Hou/PRO/PRO_co_quantified_sample.csv")
write.csv(RNA_co_quantified_sample,file = "./Hou/RNA/RNA_co_quantified_sample.csv")

rm(RNA_sample)
rm(PRO_sample)
rm(temp)
#log#
PRO_co_quantified_sample <- as.matrix(PRO_co_quantified_sample)
PRO_co_quantified_sample_log = preproc.filter.fn(PRO_co_quantified_sample,dat.type = "intensity", filter.col.prp=1)	#进行对数的计算以及过滤
#log#
RNA_co_quantified_sample <- as.matrix(RNA_co_quantified_sample)
RNA_co_quantified_sample_log = preproc.filter.fn(RNA_co_quantified_sample,dat.type = "TPM or RPKM", proc.zero ="perturbed by a small value",exp.thres=0)

Sample_info <- read.csv("./Hou/Sample_info_Sample.csv")
PRO_tiss.abd <- tiss.abd.fn(PRO_co_quantified_sample_log, Sample_info)
RNA_tiss.abd <- tiss.abd.fn(RNA_co_quantified_sample_log, Sample_info)
PRO_tiss.abd <- PRO_tiss.abd[, c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")]
RNA_tiss.abd <- RNA_tiss.abd[, c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")]

#calculate TS scores #
PRO_result = AdaTiSS(PRO_co_quantified_sample_log,tiss.abd = PRO_tiss.abd)
RNA_result = AdaTiSS(RNA_co_quantified_sample_log,tiss.abd = RNA_tiss.abd,adjust.opt = 0)
#get TS score#
PRO_Score_Sample <- PRO_result$ada.s
RNA_Score_Sample <- RNA_result$ada.s
PRO_Score_Tissue <- PRO_result$ada.z
RNA_Score_Tissue_raw <- RNA_result$ada.z
RNA_Score_Tissue <- RNA_result$ada.z

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
####TSscore filter ####
##RNA filter 1:

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

row.names(RNA_raw_tiss_abd) <- row.names(RNA_co_quantified_sample)
colnames(RNA_raw_tiss_abd) <- col_names
row.names(PRO_raw_tiss_abd) <- row.names(PRO_co_quantified_sample)
colnames(PRO_raw_tiss_abd) <- col_names

RNA_raw_tiss_abd <- apply(RNA_raw_tiss_abd, 2, function(x) ifelse(is.na(x), 0, x))
RNA_filter_1 <- RNA_raw_tiss_abd[rowSums(RNA_raw_tiss_abd< 1) == ncol(RNA_raw_tiss_abd), ]
LIST <- rownames(RNA_filter_1)
RNA_Score_Tissue[LIST,] <- "NA_all_tissues_tpm_less_1"
rm(RNA_filter_1)
##RNA filter2

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
####classify####
####protein level classify ####
PRO_filtered <- PRO_Score_Tissue
#PRO_filtered[grepl("NA", PRO_filtered)] <- NA
PRO_filtered <- as.data.frame(PRO_filtered)
##housekeeping
PRO_housekeeping_raw <- PRO_raw_tiss_abd[rowSums(is.na(PRO_raw_tiss_abd)) == 0, ]
housekeeping_raw_list<- row.names(PRO_housekeeping_raw)
PRO_housekeeping_raw <-subset(PRO_filtered, rownames(PRO_filtered) %in% housekeeping_raw_list) 
PRO_housekeeping_raw <- as.data.frame(PRO_housekeeping_raw)
PRO_housekeeping <-PRO_housekeeping_raw[rowSums(is.na(PRO_housekeeping_raw)) == 0 & apply(PRO_housekeeping_raw, 1, function(row) all(row < 2)), ]

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
####RNA level classify ####
RNA_filtered <- RNA_Score_Tissue

RNA_filtered[grepl("NA", RNA_filtered)] <- NA
RNA_filtered <- as.data.frame(RNA_filtered)
for (col in 1:ncol(RNA_filtered)) {

  RNA_filtered[, col] <- as.numeric(RNA_filtered[, col])
}
#RNA Housekeeping
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

write.csv(RNA_Enriched,file = "./Hou/RNA/RNA_Enriched.csv")
write.csv(RNA_housekeeping,file = "./Hou/RNA/RNA_housekeeping.csv")
write.csv(RNA_Enriched_not_specific,file = "./Hou/RNA/RNA_Enriched_not_specific.CSV")
write.csv(RNA_Specific,file = "D:/Workspace/Glycine_max/Tissue_Specifity/Hou/RNA/RNA_Specific.csv")
write.csv(RNA_others,file = "./Hou/RNA/RNA_others.csv")
