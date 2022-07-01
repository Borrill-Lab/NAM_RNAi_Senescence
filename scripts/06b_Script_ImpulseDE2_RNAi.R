# Aim is to run ImpulseDE2 on RNAi timecourse to compare to gradient tool results

# Tayyaba Andleeb 
# 9.07.2021

# install ImpulseDE2
source("https://bioconductor.org/biocLite.R")
biocLite("ImpulseDE2")

library(ImpulseDE2)

out_dir <- "C:/Users/TXA013/3D Objects/Borril Lab/ImpulseDE2/ImpulseDE2_RNAi/"

# this is the example from the tutorial using a simulated dataset
lsSimulatedData <- simulateDataSetImpulseDE2(
  vecTimePointsA   = rep(seq(1,8),3),
  vecTimePointsB   = NULL,
  vecBatchesA      = NULL,
  vecBatchesB      = NULL,
  scaNConst        = 30,
  scaNImp          = 10,
  scaNLin          = 10,
  scaNSig          = 10,
  scaMuBatchEffect = NULL,
  scaSDBatchEffect = NULL,
  dirOutSimulation = NULL)

# the count data is in the form:
head(lsSimulatedData$matObservedCounts)

# the annotation is in the form:
head(lsSimulatedData$dfAnnotation)

objectImpulseDE2 <- runImpulseDE2(
  matCountData    = lsSimulatedData$matObservedCounts, 
  dfAnnotation    = lsSimulatedData$dfAnnotation,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 1 )

head(objectImpulseDE2$dfImpulseDE2Results)


### now want to do it with my RNAi dataset
# should use counts as input
# only want to use genes >0.5 tpm

## do it for RNAi first
count.data <- read.csv(file="C:/Users/TXA013/3D Objects/Borril Lab/Source_Data/RNAi_timecourse_count.tsv", sep = "\t")
head(count.data)
dim(count.data)

# just use FLB data for 7 timepoints
counts <- count.data[,c(1:21)]
head(counts)
print("dimensions of count data before filtering")
dim(counts)
colnames(counts)

# now select only genes with > 0.5 tpm in at least 1 timepoint FLB
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

# clean up workspace to remove unnecessary dataframes
rm(tpmData_FLB_7timepoints)

# merge together counts_conf  and tpmData_FLB_av
counts_max_FLB <- merge(counts, tpmData_FLB_av, by.x = 0, by.y = 0)
head(counts_max_FLB)
dim(counts_max_FLB)
colnames(counts_max_FLB)
# select only rows with a maxtpm >0.5
counts_FLB <- counts_max_FLB[which(counts_max_FLB$maxtpm>0.5),]


## remove LC genes
head(counts_FLB)
counts_FLB_7timepoint <- counts_FLB[!grepl("LC",counts_FLB$Row.names), ]
head(counts_FLB_7timepoint)
dim(counts_FLB_7timepoint)

# make rownames correct
rownames(counts_FLB_7timepoint) <- counts_FLB_7timepoint[,1]
counts_FLB_7timepoint <- counts_FLB_7timepoint[,-1]
head(counts_FLB_7timepoint)
#head(row.names(counts_FLB))

# remove tpm and maxtpm columns
counts_FLB_7timepoint <- counts_FLB_7timepoint[,1:21]
head(counts_FLB_7timepoint)
dim(counts_FLB_7timepoint)

# re-organise
counts_FLB_ordered <- counts_FLB_7timepoint[,c(1,4,7,10,13,16,19,2,5,8,11,14,17,20,3,6,9,12,15,18,21)]
head(counts_FLB_ordered)

### now make annotation dataframe with sample details. Need 4 columns "Sample", "Condition", "Time" and "Batch")
annotationdf <- data.frame("Sample" = colnames(counts_FLB_ordered), 
                           "Condition" = rep("case",21), 
                           "Time" = rep(c(3,7,10,13,15,19,26),3),
                           "Batch" = rep("B_NULL",21))
head(annotationdf)
annotationdf

# now let's run ImpulseDE2 on just 15 genes as a test

test_counts_FLB_ordered <- as.matrix(round(counts_FLB_ordered[1:15,]))
test_counts_FLB_ordered

objectImpulseDE2 <- runImpulseDE2(
  matCountData    = test_counts_FLB_ordered, 
  dfAnnotation    = annotationdf,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 1 )

head(objectImpulseDE2$dfImpulseDE2Results)

# let's plot the gene expression
library(ggplot2)
lsgplotsGenes <- plotGenes(
  vecGeneIDs       = NULL,
  scaNTopIDs       = 3,
  objectImpulseDE2 = objectImpulseDE2,
  boolCaseCtrl     = FALSE,
  dirOut           = NULL,
  strFileName      = NULL,
  vecRefPval       = NULL, 
  strNameRefMethod = NULL,
  boolSimplePlot=FALSE)
print(lsgplotsGenes[[1]])


### now run for real on full set of ~52,626 genes
head(counts_FLB_ordered)
dim(counts_FLB_ordered)

counts_FLB_ordered <- as.matrix(round(counts_FLB_ordered))
head(counts_FLB_ordered)

real_objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts_FLB_ordered, 
  dfAnnotation    = annotationdf,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 1 )

head(real_objectImpulseDE2$dfImpulseDE2Results)

# save output
write.csv(file=paste0(out_dir,"RNAi_ImpulseDE2_results.csv"), real_objectImpulseDE2$dfImpulseDE2Results, row.names = F)


result_df <- real_objectImpulseDE2$dfImpulseDE2Results
head(result_df)

head(result_df[result_df$padj <0.05,])
nrow(result_df[result_df$padj <0.05,])

nrow(result_df[result_df$padj <0.01,])

nrow(result_df[result_df$padj <0.001,])

## for some reason making plots doesn't show the dots of the normalised counts ###
library(ggplot2)
lsgplotsGenes <- plotGenes(
  vecGeneIDs       = NULL,
  scaNTopIDs       = 5,
  objectImpulseDE2 = real_objectImpulseDE2,
  boolCaseCtrl     = FALSE,
  dirOut           = NULL,
  strFileName      = NULL,
  vecRefPval       = NULL, 
  strNameRefMethod = NULL,
  boolSimplePlot=FALSE)

print(lsgplotsGenes[[1]]) 


# let's compare ImpulseDE2 to gradient tool results
# get gradient tool results
gradienttool_res <- read.csv(file="C:/Users/TXA013/3D Objects/Borril Lab/Gradient Tool/RNAi_Gradient Tool/sorted_up_down_data_HC_only_RNAi.csv")
head(gradienttool_res)
dim(gradienttool_res)

gradienttool_res[gradienttool_res$X =="TraesCS4B02G311100",]

# now select only DE genes from gradient tool
gradienttool_res_DE <- gradienttool_res[gradienttool_res$pattern != "0_0_0_0_0_0_0",]
head(gradienttool_res_DE)
dim(gradienttool_res_DE)



## now try padj <0.01
# now select only DE genes from ImpulseDE2
impulse_DE <- result_df[result_df$padj <0.01,c(1,3)]
head(impulse_DE)
dim(impulse_DE)

# merge together and find overlap
merged_res <- merge(gradienttool_res_DE, impulse_DE, by.x="X", by.y="Gene", all = T)
head(merged_res)

dim(merged_res)
sum(complete.cases(merged_res))

## now try padj <0.05
# now select only DE genes from ImpulseDE2
impulse_DE <- result_df[result_df$padj <0.05,c(1,3)]
head(impulse_DE)
dim(impulse_DE)

# merge together and find overlap
merged_res <- merge(gradienttool_res_DE, impulse_DE, by.x="X", by.y="Gene", all = T)
head(merged_res)

dim(merged_res)
sum(complete.cases(merged_res))

# now select only DE genes from ImpulseDE2 padj0.001
impulse_DE <- result_df[result_df$padj <0.001,c(1,3)]
head(impulse_DE)
dim(impulse_DE)

# merge together and find overlap
merged_res <- merge(gradienttool_res_DE, impulse_DE, by.x="X", by.y="Gene", all = T)
head(merged_res)

dim(merged_res)
sum(complete.cases(merged_res))

# are those which are not shared the ones with lower padj?
head(merged_res)

write.csv(file=paste0(out_dir,"merged_impuseDE2_padj0.001_gradient_tool_RNAi.csv"),merged_res,row.names = F)
