library(MASS)
options(stringsAsFactors = F)
library(relaimpo)
library(tidyverse)
library(ggplot2)
library(ggcor)
library(ggthemes)
PRO <- read.csv("D:/Workspace/Glycine_max/m6A/Distribution_Box/PRO_Raw_Gene_Tissue.csv")
Tissue_Info <- colnames(PRO[,c(2:15)])
TPM <- read.csv("D:/Workspace/Glycine_max/m6A/Distribution_Box/RNA_Tissue_1.csv",row.names = 1)
TPM <- log2(TPM)
TPM$ID <- rownames(TPM)
TPM <- TPM[,c(15,seq(1:14))]
Intensity <- read.csv("D:/Workspace/Glycine_max/m6A/3_19/m6A.logOR.csv/m6A_Intensity.csv",row.names = 1)
Intensity$Peak.ID <- rownames(Intensity)
Intensity <- Intensity[,c(15,seq(1:14))]
Intensity$Peak.ID <- sub("^[^.]*\\.", "", Intensity$Peak.ID)
C_m6A <- read.delim("D:/Workspace/Glycine_max/m6A/3_19/conserved.txt")
C_m6A <- C_m6A[,c(6,13)]
colnames(C_m6A) <- c("Peak.ID","ID")
m6A <- C_m6A
folder_path <-"D:/Workspace/Glycine_max/m6A/3_19/nonconserved/nonconserved/"
for (i in 1:14)
{
  
  file_path <- file.path(folder_path, paste0("non", i, "conserved.anno.peak.new.anno.txt"))
  # load data
  Data <- read.delim(file_path)
  Data <- Data[,c(6,14)]
  colnames(Data)<- colnames(C_m6A)
  m6A <- rbind(m6A,Data)
}
m6A <- subset(m6A, !duplicated(m6A))
m6A_Data <- left_join(Intensity,m6A,by="Peak.ID")
m6A_Data <- m6A_Data[!is.na(m6A_Data$ID),]
m6A_Data <- m6A_Data[,-1]
m6A_Data <- m6A_Data[,c(15,seq(1:14))]
rownames(m6A_Data) <- seq(1:nrow(m6A_Data))
library(dplyr)
m6A_Data <- m6A_Data %>%
  group_by(ID) %>%
  summarize(across(.cols = ISD:MSD_R8, .fns = sum))
rm(m6A,Data,C_m6A)
rm(Intensity)
ALL_Info <-read.csv("../ALL_Info.csv")

for (i in 5:25) {
  ALL_Info [, i] <- ALL_Info [, i]*3 / ALL_Info$CDS_length
}
PRO_Isform <- read.csv("D:/Workspace/Glycine_max/m6A/Distribution_Box/PRO_Raw_Isform_Tissue.csv")
ALL_Info <- subset(ALL_Info,ALL_Info$ID %in% PRO_Isform$X)
ALL_Info$ID <- gsub("\\..*", "", ALL_Info$ID)
ALL_Info <- ALL_Info[!duplicated(ALL_Info$ID), ]
ALL_Info <- ALL_Info[,-4]
ALL_Info$Gene_length <- as.numeric(ALL_Info$Gene_length)
List <- intersect(intersect(PRO$ID,TPM$ID),m6A_Data$ID)
PRO <- subset(PRO,PRO$ID %in% List)
TPM<-subset(TPM,TPM$ID %in% List)
m6A_Data<-subset(m6A_Data,m6A_Data$ID %in% List)
rm(PRO_Isform)
for (i in 2:15) 
  {
    Name_PRO <- PRO[!is.na(PRO[,i]),1]
    Name_RNA <- TPM[!is.na(TPM[,i]),1]
    Name_m6A <- m6A_Data[m6A_Data[,i]>0,1]
    Name<-intersect(Name_PRO,Name_RNA)
    PRO_Temp <- subset(PRO,PRO$ID %in% Name)
    Data1<- as.data.frame(PRO_Temp[,i])
    colnames(Data1) <- "PRO"
    RNA_Temp <- subset(TPM,TPM$ID %in% Name)
    Data2<- as.data.frame(RNA_Temp[,i])
    colnames(Data2) <- "RNA_abundance"
    m6A_Temp <- subset(m6A_Data,m6A_Data$ID %in% Name)
    Data3<- m6A_Temp[,i]
    colnames(Data3) <- "m6A"
    Info_Temp <-subset(ALL_Info,ALL_Info$ID %in% Name)
    Data <- cbind(Data1,Data3,Data2,Info_Temp)
    rownames(Data) <- Name
    library(MASS)
    Data <- Data[,-4]
    Test <- Data[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26)]
    full.model <- lm(PRO ~., data =Test)
    #step.model <- stepAIC(full.model, direction = "both",trace = FALSE)
    #summary(step.model)
    #relaimpo::calc.relimp.lm(full.model,rela=T)
    importance <- relaimpo::calc.relimp.lm(full.model,rela=T)
    Temp <- data.frame(factor=names(importance@lmg),importance=importance@lmg)
    assign(Tissue_Info[i-1], Temp)
    } 


for (i in 1:14) 
{
  name <- Tissue_Info[i]
  Data <- get(name)
  colnames(Data)<- c("factor",Tissue_Info[i])
  name <- Tissue_Info[i]
  assign(name, Data)
}
importance_Data <- cbind(ISD$factor,ISD$ISD,CT_VE$CT_VE,RT_VE$RT_VE,
                         ULF_V1$ULF_V1,SM_V1$SM_V1,RT_V1$RT_V1,CT_V1$CT_V1,
                         TLF_V1$TLF_V1,FL_R2$FL_R2,RT_R5$RT_R5,Pod_R6$Pod_R6,
                         RTN_R5$RTN_R5,GSD_R7$GSD_R7,MSD_R8$MSD_R8)
colnames(importance_Data) <- c("Factor",Tissue_Info)

importance_Data <- as.data.frame(importance_Data)
for (i in 2:ncol(importance_Data)) {
  importance_Data[,i] <- as.numeric(importance_Data[,i])
}
rownames(importance_Data) <- importance_Data$Factor
importance_Data <- importance_Data[,2:15]

Codon_usage <- colSums(importance_Data[5:24,])
importance_Data <- rbind(importance_Data[1:4,],Codon_usage)
row.names(importance_Data) <- c("m6A","RNA abundance","Gene length","Exon number","Codon usage")


Data <- data.frame(matrix("", nrow = 70, ncol =3))
colnames(Data) <- c("value","Tissue","factor")
Data$Tissue <- rep(Tissue_Info,time=5)
Data$factor <- rep(rownames(importance_Data),each=14)
Data[1:14,1] <- as.numeric(importance_Data[1,])
Data[15:28,1] <- as.numeric(importance_Data[2,])
Data[29:42,1] <- as.numeric(importance_Data[3,])
Data[43:56,1] <- as.numeric(importance_Data[4,])
Data[57:70,1] <- as.numeric(importance_Data[5,])
Data$value <- as.numeric(Data$value)
Data$factor <- factor(Data$factor,levels = rev(c("RNA abundance","Codon usage", "m6A", "Gene length","Exon number")))


write.csv(Data,file="../Data.csv")
####violin plot ####
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(cowplot)
library(reshape2)
library(tidyverse)
library(ggsci)
ggplot(Data,aes(factor,value))+
  geom_boxplot(notch = F,aes(fill=factor),color="black",outlier.shape = NA)+coord_flip()+
  stat_boxplot(geom = "errorbar",aes(ymin = ..ymax..,color=factor), width = 0.3, size = .5) +
  stat_boxplot(geom = "errorbar",aes(ymax = ..ymin..,color=factor), width = 0.3, size = .5) +
  scale_fill_manual(values = rev(c("#EF3B2C","#4292C6","#41AB5D","#807DBA","#F16913")))+
  scale_color_manual(values = rev(c("#EF3B2C","#4292C6","#41AB5D","#807DBA","#F16913")))+
  theme(legend.title = element_text(colour = NA))+
  theme_test() + theme(legend.position = "none")+
  labs(y="Relative contribution to protein abundance",x="")
 
