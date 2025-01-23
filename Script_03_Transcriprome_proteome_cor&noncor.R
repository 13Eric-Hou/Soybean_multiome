####Core-RNA correlation####
rm(list = ls())
RNA <- read.csv("D:/Workspace/Glycine_max/RNA_Seq/TPM/Gene_TPM.csv",row.names = 1)
RNA_Tissue <- data.frame(matrix("", nrow = 56750, ncol =14))
rownames(RNA_Tissue) <- rownames(RNA)
Tissue_info <-  c("ISD","CT_VE","RT_VE","ULF_V1","SM_V1","RT_V1","CT_V1","TLF_V1","FL_R2","RT_R5","Pod_R6","RTN_R5","GSD_R7","MSD_R8")
colnames(RNA_Tissue) <- Tissue_info
for (i in 1:nrow(RNA))
  {
    for (j in 1:14) 
    {
      RNA_Tissue[i,j] <- mean(as.numeric(RNA[i,c((3*j-2),(3*j-1),(3*j))]),na.rm = T)  
    }
   
}
for (i in 1:14) 
  {
    RNA_Tissue[,i] <- as.numeric(RNA_Tissue[,i])
  }
RNA_Tissue[is.nan(as.matrix(RNA_Tissue))] <- NA
##adjust core and non_core
Temp <- RNA_Tissue
RNA_Tissue <- Temp
RNA_Tissue[RNA_Tissue >= 1] <- "Core"
RNA_Tissue[RNA_Tissue < 1] <- "Non-Core"
RNA_Expression <-Temp
RNA_Expression<- log2(RNA_Expression)

#RNA_Expression[is.na(RNA_Expression)] <- 0

##correlation analysis（Pairwise global Pearson’s correlation）##
Result_1 <- data.frame(matrix("", nrow = 14, ncol =14))
colnames(Result_1) <- Tissue_info
rownames(Result_1) <- Tissue_info
for (i in 1:14) 
  {
    for (j in 1:14) 
    {
      ID_1<- rownames(RNA_Tissue[RNA_Tissue[,i]=="Core",])
      ID_2<- rownames(RNA_Tissue[RNA_Tissue[,j]=="Core",])
      J <- intersect(ID_1,ID_2)
      Expression_1 <- subset(RNA_Expression,rownames(RNA_Expression) %in% J)
      Expression_1 <- Expression_1[!is.na(Expression_1[,i]),]
      Expression_2 <- subset(RNA_Expression,rownames(RNA_Expression) %in% J)
      Expression_2 <- subset(Expression_2,rownames(Expression_2) %in% rownames(Expression_1))

      Expression_2 <- Expression_2[!is.na(Expression_2[,j]),]
      Expression_1 <- subset(Expression_1,rownames(Expression_1) %in% rownames(Expression_2))
      Result_1[i,j] <- cor(Expression_1[,i],Expression_2[,j],method = "pearson",use = "pairwise.complete.obs")
    }
   
  }
##Non-core correlation analysis（Pairwise global Pearson’s correlation）##
Result_2 <- data.frame(matrix("", nrow = 14, ncol =14))
colnames(Result_2) <- Tissue_info
rownames(Result_2) <- Tissue_info
for (i in 1:14) 
{
  for (j in 1:14) 
  {
    ID_1<- rownames(RNA_Tissue[RNA_Tissue[,i]=="Non-Core",])
    ID_2<- rownames(RNA_Tissue[RNA_Tissue[,j]=="Non-Core",])
    J <- intersect(ID_1,ID_2)
    Expression_1 <- subset(RNA_Expression,rownames(RNA_Expression) %in% J)
    Expression_2 <- subset(RNA_Expression,rownames(RNA_Expression) %in% J)
    Expression_1 <- Expression_1[!is.na(Expression_1[,i]),]
    Expression_2 <- subset(Expression_2,rownames(Expression_2) %in% rownames(Expression_1))
    Expression_2 <- Expression_2[!is.na(Expression_2[,j]),]
    Expression_1 <- subset(Expression_1,rownames(Expression_1) %in% rownames(Expression_2))
    Result_2[i,j] <- cor(Expression_1[,i],Expression_2[,j],method = "pearson",use = "pairwise.complete.obs")
  }
  
}
row_order <- c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
               "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2")
col_order <- c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
               "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2")
Result_1 <- Result_1[row_order,]
Result_1 <-Result_1[,col_order]
Result_2 <- Result_2[row_order,]
Result_2 <-Result_2[,col_order]
Result <- Result_1
for (i in 2:14) 
  {
    Result[i:14,i-1] <- Result_2[i:14,i-1]
}
##pheatmap##
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

for (i in 1:14) 
  {
    Result[,i] <-as.numeric(Result[,i] )
}
Result <- as.matrix(Result)

pheatmap(Result,
         fontsize = 15,
         cellwidth =15,cellheight = 15,
         cluster_cols = F,cluster_rows = F,
         breaks = c(seq(0, 0.4, length.out = 50), 
                    seq(0.401,1, length.out = 50)))

####protein expression####
PRO_Expression <-read.csv("D:/Workspace/Glycine_max/Data/PRO_VSN_Gene_Tissue.csv",row.names = 1)
Result_1 <- data.frame(matrix("", nrow = 14, ncol =14))
colnames(Result_1) <- Tissue_info
rownames(Result_1) <- Tissue_info
for (i in 1:14) 
{
  for (j in 1:14) 
  {
    ID_1<- rownames(RNA_Tissue[RNA_Tissue[,i]=="Core",])
    ID_2<- rownames(RNA_Tissue[RNA_Tissue[,j]=="Core",])
    J <- intersect(ID_1,ID_2)
    Expression_1 <- subset(PRO_Expression,rownames(PRO_Expression) %in% J)
    Expression_2 <- subset(PRO_Expression,rownames(PRO_Expression) %in% J)
    Expression_1 <- Expression_1[!is.na(Expression_1[,i]),]
    Expression_2 <- subset(Expression_2,rownames(Expression_2) %in% rownames(Expression_1))
    Expression_2 <- Expression_2[!is.na(Expression_2[,j]),]
    Expression_1 <- subset(Expression_1,rownames(Expression_1) %in% rownames(Expression_2))
    Result_1[i,j] <- cor(Expression_1[,i],Expression_2[,j],method = "pearson",use = "pairwise.complete.obs")
  }
  
}
##Non-core gene correlation analysis（Pairwise global Pearson’s correlation）##
Result_2 <- data.frame(matrix("", nrow = 14, ncol =14))
colnames(Result_2) <- Tissue_info
rownames(Result_2) <- Tissue_info
for (i in 1:14) 
{
  for (j in 1:14) 
  {
    ID_1<- rownames(RNA_Tissue[RNA_Tissue[,i]=="Non-Core",])
    ID_2<- rownames(RNA_Tissue[RNA_Tissue[,j]=="Non-Core",])
    J <- intersect(ID_1,ID_2)
    Expression_1 <- subset(PRO_Expression,rownames(PRO_Expression) %in% J)
    Expression_2 <- subset(PRO_Expression,rownames(PRO_Expression) %in% J)
    Expression_1 <- Expression_1[!is.na(Expression_1[,i]),]
    Expression_2 <- subset(Expression_2,rownames(Expression_2) %in% rownames(Expression_1))
    Expression_2 <- Expression_2[!is.na(Expression_2[,j]),]
    Expression_1 <- subset(Expression_1,rownames(Expression_1) %in% rownames(Expression_2))
    Result_2[i,j] <- cor(Expression_1[,i],Expression_2[,j],method = "pearson",use = "pairwise.complete.obs")
  }
  
}
row_order <- c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
               "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2")
col_order <- c("Pod_R6","GSD_R7","MSD_R8","ISD","RT_VE","RT_V1","RT_R5","RTN_R5",
               "CT_VE","CT_V1","ULF_V1","TLF_V1","SM_V1","FL_R2")
Result_1 <- Result_1[row_order,]
Result_1 <-Result_1[,col_order]
Result_2 <- Result_2[row_order,]
Result_2 <-Result_2[,col_order]
Result <- Result_1
for (i in 2:14) 
{
  Result[i:14,i-1] <- Result_2[i:14,i-1]
}
##pheatmap##
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

for (i in 1:14) 
{
  Result[,i] <-as.numeric(Result[,i] )
}
Result <- as.matrix(Result)
coul <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

pheatmap(Result,
         fontsize = 15,
         cellwidth =15,cellheight = 15,
         cluster_cols = F,cluster_rows = F,
         breaks = c(seq(0, 0.4, length.out = 50), 
                    seq(0.401, 0.9, length.out = 50)))
