rm(list=ls())
library(GSVA)
library(GSEABase)
setwd("C:/Users/penguin/Desktop/学习资料/功能基因组/lab2-免疫细胞浸润分析/lab2-免疫细胞浸润分析/homework")

#library(CIBERSORT)
#data(LM22)
#dim(LM22)
#LM22[1:2,]

library(readxl)
lm22_data<-read_excel("41592_2015_BFnmeth3337_MOESM207_ESM.xls",sheet = 4)
dim(lm22_data)
lm22_matrix <- as.matrix(lm22_data[2:549,])
gene <- lm22_matrix[2:548,1]
cell <- lm22_matrix[1,2:23]

lm22_matrix <- lm22_matrix[2:548,2:23]
colnames(lm22_matrix) <- cell
rownames(lm22_matrix) <- gene

lm22 <-c()
for(i in 1:dim(lm22_matrix)[2]){
  genes <- rownames(lm22_matrix)[which(lm22_matrix[,i] == 1)]
  length(genes)<-70
  genes <- c(NA,genes)
  lm22 <- cbind(lm22,genes)
}
lm22 <- t(lm22)
write.table(x=lm22,file ="LM22.gmt",quote = FALSE, sep = "\t",col.names = F, row.names = cell)
save.image("LM22.Rdata")


#TCGA
library(BiocManager)
library(TCGAbiolinks)
library(stringr)
library(tidyverse)

CancerProject <- "TCGA-GBM"
DataDirectory <- paste0("./",gsub("-","_",CancerProject)) 
# gsub("-","_",CancerProject) # 替换"-"为"_"

# 数据名称
FileNameData <- paste0(DataDirectory, "_","expression",".rda")

query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Transcriptome Profiling",
                  experimental.strategy = "RNA-Seq",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")

samplesDown <- getResults(query,cols=c("cases"))
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP") #C

dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT") #正常实体组织样品编号
dataSmTR <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TR") #复发实体组织样品编号

dataSmTP <- gsub("-",".",dataSmTP)
dataSmNT <- gsub("-",".",dataSmNT)
dataSmTR <- gsub("-",".",dataSmTR)

GDCdownload(query)
data <- GDCprepare(query = query, 
                   save = TRUE, 
                   save.filename = FileNameData)

#准备表达矩阵
GBMRnaseqSE <- GDCprepare(query)
library(SummarizedExperiment)
assayNames(data)

rowdata <- rowData(data)  # 查看rowdata内容
names(rowdata)
table(rowdata$gene_type) # 统计基因类型，下面可以提取的选择就在这里

head(rowdata$gene_name,10)   # gene_name就是symbol
length(rowdata$gene_name) #有多少gene
data_miRNA <- data[rowdata$gene_type == "miRNA"]  # 提取miRNA
data_mRNA <- data[rowdata$gene_type == "protein_coding",] # 提取mRNA


#####三、提取表达矩阵#####

# 提取表达矩阵
# mRNA的counts矩阵
test.mRNA.counts <- assay(data_mRNA,"unstranded") # counts矩阵用于后续差异分析

# mRNA的tpm矩阵
test.mRNA.tpm <- assay(data_mRNA,"tpm_unstrand")

# mRNA的fpkm矩阵
test.mRNA.fpkm <- assay(data_mRNA,"fpkm_unstrand")

# 添加gene_symbol,基因名称
mRNA.symbol <- rowData(data_mRNA)$gene_name

# 合并
test.mrna.frame <- cbind(as.data.frame(mRNA.symbol),
                         as.data.frame(test.mRNA.counts))

dim(test.mrna.frame) #基因/样品数量

#####四、初步处理#####

# 去重数据
qc = as.matrix(test.mrna.frame)
rownames(qc)=qc[,1] # gene_symbol
exp=qc[,2:ncol(qc)] # matrix
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

# 对重复的基因取平均值
#BiocManager::install("limma")
library(limma)
data = limma::avereps(data) # 去重,出现多行取平均值

# 过滤表达量低的基因,可根据基因数目进行调整
CHOL.test=data[rowMeans(data)>5,]
output <- rbind(colnames(CHOL.test),CHOL.test)
output <- output[-1,]
dim(output)
write.csv(output, file="GBM_mRNA矩阵去重.csv",quote = FALSE)
write.table(output,file ="GBM_mrna_expression.gmt",quote = FALSE, sep = "\t",col.names = T, row.names = T)
save.image("TCGA_data_clean.Rdata")


#GSVA
library(GSEABase)
library(GSVA)
geneSet <- getGmt("./LM22.gmt") 
class(geneSet)

gene_exp = read.table("./GBM_mrna_expression.gmt", sep="\t", header=TRUE, row.names=1) #读入表达谱
gene_exp = as.matrix(gene_exp)

gsva_matrix<- gsva(gene_exp, geneSet, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

gsva_matrix[1:2,]
write.csv(gsva_matrix,"./GBM_ssGSEA_score.csv")

#4.画热图
library(pheatmap)
gsva_matrix1<- t(scale(t(gsva_matrix)))#默认按照列进行归一化
# gsva_matrix1[gsva_matrix1< -2] <- -2
# gsva_matrix1[gsva_matrix1>2] <- 2



normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}#设定normalization函数,将值映射到[0，1]之间

nor_gsva_matrix1 <- normalization(gsva_matrix1)
annotation_col <- c(rep("Other",ncol(gene_exp)))
for(i in 1: ncol(gene_exp)){
  if(colnames(gene_exp)[i] %in% dataSmTP){
    annotation_col[i] = "Primary Solid Tumor"
  }
  else if(colnames(gene_exp)[i] %in% dataSmNT){
    annotation_col[i] = "Solid Tissue Normal"
  }
  else if(colnames(gene_exp)[i] %in% dataSmTR){
    annotation_col[i] = "Recurrent Solid Tumor"
  }
}

annotation_col <- data.frame(phenotype = annotation_col)
rownames(annotation_col)<-colnames(gene_exp)#使编号能互相对应
bk = unique(c(seq(0,1, length=100)))#设定热图参数
pheatmap(nor_gsva_matrix1,
         color=colorRampPalette(c("blue","white","deeppink"))(100),
         show_colnames = F,
         cluster_rows = T,cluster_cols = T,
         annotation_col = annotation_col,
         breaks=bk,cellwidth=5,cellheight=8,
         fontsize=5,
         filename = './GBM_ssgsea.pdf',width = 16)#画热图


#生存分析
#TCGA clinic 数据处理
GBM_clinic <- GDCquery_clinic(project = "TCGA-GBM", type = "clinical")
dim(GBM_clinic)
#nor_gsva_matrix1富集矩阵

newid<-lapply(strsplit(colnames(data),'-'),function(i){paste0(i[1:3],collapse = '-')})
newid<-sapply(1:ncol(data),function(i){newid[[i]]})

nor_gsva_matrix1<-rbind(nor_gsva_matrix1,newid)

colnames(GBM_clinic)[2]<-'newid'
df_OS<-merge(GBM_clinic,t(nor_gsva_matrix1),by="newid")


alive <- df_OS[which(df_OS$vital_status=="Alive"),]
dead <- df_OS[which(df_OS$vital_status=="Dead"),]
alive$time <- alive$days_to_last_follow_up
dead$time <- dead$days_to_death

all_patients <- rbind(alive,dead)

patient_survive <- cbind(all_patients$newid,all_patients$time,all_patients$vital_status)
colnames(patient_survive) <- c("newid","time","status")

GBM_survive_matrix <- merge(t(nor_gsva_matrix1),patient_survive,by="newid")

library(limma)
GBM_survive_matrix= limma::avereps(GBM_survive_matrix,ID=GBM_survive_matrix$newid)
GBM_survive_matrix[GBM_survive_matrix == "Alive"] = 0
GBM_survive_matrix[GBM_survive_matrix == "Dead"] = 1
rownames(GBM_survive_matrix) <- GBM_survive_matrix[,1]
GBM_survive_matrix<-GBM_survive_matrix[,-1]
GBM <- as.data.frame(GBM_survive_matrix)

#生存分析
library(survival)
library(survminer)

B_cell_naive <- c()
for(i in 1: nrow(GBM)){
  if(GBM[i,1] >= mean(as.numeric(GBM[,1]))){
    B_cell_naive[i] = "high"
  }else B_cell_naive[i] = "low"
}
GBM <- cbind(GBM,B_cell_naive)
fit1 <- survfit(Surv(as.numeric(time),as.numeric(status)) ~ B_cell_naive , data = GBM)

Macrophages_M2 <- c()
for(i in 1: nrow(GBM)){
  if(GBM[i,1] >= mean(as.numeric(GBM[,16]))){
    Macrophages_M2[i] = "high"
  }else Macrophages_M2[i] = "low"
}
GBM <- cbind(GBM,Macrophages_M2)
fit2 <- survfit(Surv(as.numeric(time),as.numeric(status)) ~ Macrophages_M2 , data = GBM)

ggsurvplot(
  fit2,                     
  data = GBM,             
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = TRUE,         
  xlim = c(0,500),         # 横坐标轴范围，相当于局部放大
  xlab = "Time in days",   # 横坐标标题
  break.time.by = 100,     # 横坐标刻度
  ggtheme = theme_light(), 
  risk.table.y.text.col = T, # risk table文字注释颜色
  risk.table.y.text = FALSE # risk table显示条形而不是文字
  )