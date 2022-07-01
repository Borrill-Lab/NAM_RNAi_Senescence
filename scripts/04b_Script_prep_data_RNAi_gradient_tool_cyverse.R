# Script to prepare data for gradient tool on cyverse to look at genes differentially expressed across time in RNA-seq timecourse
# Tayyaba Andleeb
# 08.06.2021


## read in data and only select genes with > 0.5 tpm ##

out_dir <- "C:/Users/TXA013/3D Objects/Borril Lab/Gradient Tool/RNAi_Gradient Tool/"


## do it for RNAi
## read in data and only select genes with > 0.5 tpm ##
RNAi.tpm.data <- read.csv(file="C:/Users/TXA013/3D Objects/Borril Lab/Source_Data/RNAi_timecourse_tpm.tsv", sep = "\t")
head(RNAi.tpm.data)

# select genes >0.5 tpm RNAi:
# just use FLB data
RNAi.tpm.data <- RNAi.tpm.data[,1:21]
head(RNAi.tpm.data)

# average per timepoint
RNAi.tpm.data$T3 <- (RNAi.tpm.data[,1] + RNAi.tpm.data[,2] + RNAi.tpm.data[,3]) / 3
RNAi.tpm.data$T7 <- (RNAi.tpm.data[,4] + RNAi.tpm.data[,5] + RNAi.tpm.data[,6]) / 3
RNAi.tpm.data$T10 <- (RNAi.tpm.data[,7] + RNAi.tpm.data[,8] + RNAi.tpm.data[,9]) / 3
RNAi.tpm.data$T13 <- (RNAi.tpm.data[,10] + RNAi.tpm.data[,11] + RNAi.tpm.data[,12]) / 3
RNAi.tpm.data$T15 <- (RNAi.tpm.data[,13] + RNAi.tpm.data[,14] + RNAi.tpm.data[,15]) / 3
RNAi.tpm.data$T19 <- (RNAi.tpm.data[,16] + RNAi.tpm.data[,17] + RNAi.tpm.data[,18]) / 3
RNAi.tpm.data$T26 <- (RNAi.tpm.data[,19] + RNAi.tpm.data[,20] + RNAi.tpm.data[,21]) / 3

colnames(RNAi.tpm.data[,22:28])
RNAi.tpm.data$maxtpm <- apply(RNAi.tpm.data[,22:28],1,max)
head(RNAi.tpm.data)

# now select only genes >0.5 max tpm for RNAi
RNAi_genes_0.5tpm <- RNAi.tpm.data[RNAi.tpm.data$maxtpm > 0.5,]
dim(RNAi.tpm.data)
dim(RNAi_genes_0.5tpm)
head(RNAi_genes_0.5tpm)
RNAi_genes_0.5tpm$gene <- rownames(RNAi_genes_0.5tpm)
head(RNAi_genes_0.5tpm)

#remove LC genes
RNAi_genes_0.5tpm <- RNAi_genes_0.5tpm[!grepl("LC",RNAi_genes_0.5tpm$gene), ]
head(RNAi_genes_0.5tpm)
dim(RNAi_genes_0.5tpm)


## now actually make data into format for gradient tool
reordered.FLB.RNAi <- RNAi_genes_0.5tpm[,c(1,4,7,10,13,16,19,2,5,8,11,14,17,20,3,6,9,12,15,18,21)]
head(reordered.FLB.RNAi)
dim(reordered.FLB.RNAi)

colnames(reordered.FLB.RNAi) <- c(rep(c("3", "7", "10", "13", "15", "19", "26"),3))
head(reordered.FLB.RNAi)

write.csv(file=paste0(out_dir_RNAi,"RNAi_tpm_FLB_0.5tpm_input_gradient_tool.csv"),reordered.FLB.RNAi)

