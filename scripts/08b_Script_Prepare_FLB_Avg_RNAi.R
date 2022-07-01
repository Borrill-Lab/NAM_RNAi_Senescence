# Aim is to prepare data containing Avg of TPM value and MaxTPM to identify more N responsive genes for 7 timepoints in RNAi

# Tayyaba Andleeb 
# 16.01.2022
out_dir <- "C:/Users/TXA013/3D Objects/Borril Lab/ImpulseDE2/ImpulseDE2_RNAi/"
tpmData_FLB <- read.csv(file="C:/Users/TXA013/3D Objects/Borril Lab/Source_Data/RNAi_timecourse_tpm.tsv", header=T, sep = "\t")
head(tpmData_FLB)
colnames(tpmData_FLB)
dim(tpmData_FLB)

# just use FLB data 
tpmData_FLB_7timepoints <- tpmData_FLB[,c(1:21)]
head(tpmData_FLB_7timepoints)
colnames(tpmData_FLB_7timepoints)
dim(tpmData_FLB_7timepoints)

# average per timepoint
tpmData_FLB_7timepoints$T3 <- (tpmData_FLB_7timepoints[,1] + tpmData_FLB_7timepoints[,2] + tpmData_FLB_7timepoints[,3]) / 3
tpmData_FLB_7timepoints$T7 <- (tpmData_FLB_7timepoints[,4] + tpmData_FLB_7timepoints[,5] + tpmData_FLB_7timepoints[,6]) / 3
tpmData_FLB_7timepoints$T10 <- (tpmData_FLB_7timepoints[,7] + tpmData_FLB_7timepoints[,8] + tpmData_FLB_7timepoints[,9]) / 3
tpmData_FLB_7timepoints$T13 <- (tpmData_FLB_7timepoints[,10] + tpmData_FLB_7timepoints[,11] + tpmData_FLB_7timepoints[,12]) / 3
tpmData_FLB_7timepoints$T15 <- (tpmData_FLB_7timepoints[,13] + tpmData_FLB_7timepoints[,14] + tpmData_FLB_7timepoints[,15]) / 3
tpmData_FLB_7timepoints$T19 <- (tpmData_FLB_7timepoints[,16] + tpmData_FLB_7timepoints[,17] + tpmData_FLB_7timepoints[,18]) / 3
tpmData_FLB_7timepoints$T26 <- (tpmData_FLB_7timepoints[,19] + tpmData_FLB_7timepoints[,20] + tpmData_FLB_7timepoints[,21]) / 3

# just keep the average per timepoint
#head(tpmData_FLB)
colnames(tpmData_FLB_7timepoints)[22:28]
tpmData_FLB_av <- tpmData_FLB_7timepoints[,22:28]

tpmData_FLB_av$maxtpm <- apply(tpmData_FLB_av[,1:7],1,max)
head(tpmData_FLB_av)

#  tpmData_FLB_av dataframe
TPM_max_FLB <- merge(tpmData_FLB_7timepoints[1:21], tpmData_FLB_av, by.x = 0, by.y = 0)
head(TPM_max_FLB)
dim(TPM_max_FLB)
colnames(TPM_max_FLB)
# select only rows with a maxtpm >0.5
TPM_FLB <- TPM_max_FLB[which(TPM_max_FLB$maxtpm>0.5),]


## remove LC genes
head(TPM_FLB)
TPM_FLB_7timepoint <- TPM_FLB[!grepl("LC",TPM_FLB$Row.names), ]
head(TPM_FLB_7timepoint)
dim(TPM_FLB_7timepoint)
write.csv(file=paste0(out_dir,"RNAi_FLB_Avg.csv"), TPM_FLB_7timepoint, row.names = T)