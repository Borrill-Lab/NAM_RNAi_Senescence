# Script to prepare data for gradient tool on cyverse to look at genes differentially expressed across time in RNA-seq timecourse

# Tayyaba Andleeb 
# 9.07.2021
## read in data and only select genes with > 0.5 tpm ##

out_dir <- "C:/Users/TXA013/3D Objects/Borril Lab/Gradient Tool/Control_7timepoints_Gradient Tool/"


## do it for control 
# now select only genes with > 0.5 tpm in at least 1 timepoint FLB
tpmData_FLB <- read.csv(file="C:/Users/TXA013/3D Objects/Borril Lab/Source_Data/control_timecourse_tpm.tsv", header=T, sep = "\t")
head(tpmData_FLB)
colnames(tpmData_FLB)
dim(tpmData_FLB)

# just use FLB data for 7 timepoints
tpmData_FLB <- tpmData_FLB[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,19,20,21,28,29,30)]
head(tpmData_FLB)
colnames(tpmData_FLB)
dim(tpmData_FLB)

# average per timepoint
tpmData_FLB$T3 <- (tpmData_FLB[,1] + tpmData_FLB[,2] + tpmData_FLB[,3]) / 3
tpmData_FLB$T7 <- (tpmData_FLB[,4] + tpmData_FLB[,5] + tpmData_FLB[,6]) / 3
tpmData_FLB$T10 <- (tpmData_FLB[,7] + tpmData_FLB[,8] + tpmData_FLB[,9]) / 3
tpmData_FLB$T13 <- (tpmData_FLB[,10] + tpmData_FLB[,11] + tpmData_FLB[,12]) / 3
tpmData_FLB$T15 <- (tpmData_FLB[,13] + tpmData_FLB[,14] + tpmData_FLB[,15]) / 3
tpmData_FLB$T19 <- (tpmData_FLB[,16] + tpmData_FLB[,17] + tpmData_FLB[,18]) / 3
tpmData_FLB$T26 <- (tpmData_FLB[,19] + tpmData_FLB[,20] + tpmData_FLB[,21]) / 3

# calculate the average max tpm
colnames(tpmData_FLB[,22:28])
tpmData_FLB$maxtpm <- apply(tpmData_FLB[,22:28],1,max)
head(tpmData_FLB)


# select only rows with a maxtpm >0.5
tpm_max_FLB <- tpmData_FLB[which(tpmData_FLB$maxtpm>0.5),]
head(tpm_max_FLB)
## remove LC genes
head(tpm_max_FLB)
tpm_max_FLB$Row.names <- rownames(tpm_max_FLB)
tpm_max_FLB <- tpm_max_FLB[!grepl("LC",tpm_max_FLB$Row.names), ]
head(tpm_max_FLB)
dim(tpm_max_FLB)

# remove tpm and maxtpm columns
tpm_max_FLB <- tpm_max_FLB[,1:21]
head(tpm_max_FLB)
dim(tpm_max_FLB)


# re-organise
tpm_FLB_ordered <- tpm_max_FLB[,c(1,4,7,10,13,16,19,2,5,8,11,14,17,20,3,6,9,12,15,18,21)]
head(tpm_FLB_ordered)
colnames(tpm_FLB_ordered) <- c(rep(c("3", "7", "10", "13", "15", "19", "26"),3))
head(tpm_FLB_ordered)

write.csv(file=paste0(out_dir,"control_7timepoints_tpm_FLB_0.5tpm_input_gradient_tool.csv"),tpm_FLB_ordered)

