## Run DESeq2 WT vs RNAi
# 8.2.18
# updated 31.5.2018
# updated 13.6.2022

data_dir <- "C:\\Users\\borrillp\\Documents\\Travel_May2018\\RefSeqv1.1_masked\\previous_work\\expression_per_gene\\"

setwd("C:\\Users\\borrillp\\Documents\\Travel_May2018\\RefSeqv1.1_masked\\new_May2018\\DESeq2\\")

## read in data and only select genes with > 0.5 tpm ##
RNAi.tpm.data <- read.csv(file=paste0(data_dir,"RNAi_timecourse_tpm.tsv"), sep = "\t")
head(RNAi.tpm.data)
control.tpm.data <- read.csv(file=paste0(data_dir,"control_timecourse_tpm.tsv"), sep = "\t")
head(control.tpm.data)

# select genes >0.5 tpm control:
# just use FLB data
control.tpm.data <- control.tpm.data[,1:30]
head(control.tpm.data)

# average per timepoint
control.tpm.data$T3 <- (control.tpm.data[,1] + control.tpm.data[,2] + control.tpm.data[,3]) / 3
control.tpm.data$T7 <- (control.tpm.data[,4] + control.tpm.data[,5] + control.tpm.data[,6]) / 3
control.tpm.data$T10 <- (control.tpm.data[,7] + control.tpm.data[,8] + control.tpm.data[,9]) / 3
control.tpm.data$T13 <- (control.tpm.data[,10] + control.tpm.data[,11] + control.tpm.data[,12]) / 3
control.tpm.data$T15 <- (control.tpm.data[,13] + control.tpm.data[,14] + control.tpm.data[,15]) / 3
control.tpm.data$T17 <- (control.tpm.data[,16] + control.tpm.data[,17] + control.tpm.data[,18]) / 3
control.tpm.data$T19 <- (control.tpm.data[,19] + control.tpm.data[,20] + control.tpm.data[,21]) / 3
control.tpm.data$T21 <- (control.tpm.data[,22] + control.tpm.data[,23] + control.tpm.data[,24]) / 3
control.tpm.data$T23 <- (control.tpm.data[,25] + control.tpm.data[,26] + control.tpm.data[,27]) / 3
control.tpm.data$T26 <- (control.tpm.data[,28] + control.tpm.data[,29] + control.tpm.data[,30]) / 3

# just keep the average per timepoint
colnames(control.tpm.data)[31:40]
tpmData_av_con <- control.tpm.data[,31:40]
head(tpmData_av_con)

tpmData_av_con$maxtpm <- apply(tpmData_av_con[,1:10],1,max)
head(tpmData_av_con)

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


# just keep the average per timepoint
colnames(RNAi.tpm.data)[22:28]
tpmData_av_RNAi <- RNAi.tpm.data[,22:28]
head(tpmData_av_RNAi)

tpmData_av_RNAi$maxtpm <- apply(tpmData_av_RNAi[,1:7],1,max)
head(tpmData_av_RNAi)

# now select only genes >0.5 max tpm for Control and RNAi
control_genes_0.5tpm <- tpmData_av_con[tpmData_av_con$maxtpm > 0.5,]
dim(tpmData_av_con)
dim(control_genes_0.5tpm)

RNAi_genes_0.5tpm <- tpmData_av_RNAi[tpmData_av_RNAi$maxtpm > 0.5,]
dim(tpmData_av_RNAi)
dim(RNAi_genes_0.5tpm)

merged_genes_0.5tpm <- merge(control_genes_0.5tpm, RNAi_genes_0.5tpm, by = 0, all = T)
head(merged_genes_0.5tpm)
dim(merged_genes_0.5tpm)
colnames(merged_genes_0.5tpm)[1] <- "gene"
head(merged_genes_0.5tpm)
#remove LC genes
merged_genes_0.5tpm <- merged_genes_0.5tpm[!grepl("LC",merged_genes_0.5tpm$gene), ]
head(merged_genes_0.5tpm)
dim(merged_genes_0.5tpm)

#### now read in count data and select only genes expressed >0.5 tpm ####
RNAi.count.data <- read.csv(file=paste0(data_dir,"RNAi_timecourse_count.tsv"), sep = "\t")
head(RNAi.count.data)
control.count.data <- read.csv(file=paste0(data_dir,"control_timecourse_count.tsv"), sep = "\t")
head(control.count.data)

FLB.RNAi.count.data <- RNAi.count.data[,1:21]
head(FLB.RNAi.count.data)
dim(FLB.RNAi.count.data)
FLB.RNAi.count.data.0.5tpm <- FLB.RNAi.count.data[rownames(FLB.RNAi.count.data) %in% merged_genes_0.5tpm$gene,]
head(FLB.RNAi.count.data.0.5tpm)
dim(FLB.RNAi.count.data.0.5tpm)


FLB.control.count.data <- control.count.data[,1:30]
head(FLB.control.count.data)
dim(FLB.control.count.data)
FLB.control.count.data.0.5tpm <- FLB.control.count.data[rownames(FLB.control.count.data) %in% merged_genes_0.5tpm$gene,]
head(FLB.control.count.data.0.5tpm)
dim(FLB.control.count.data.0.5tpm)

FLB.control.count.data.0.5tpm.timepoints <- FLB.control.count.data.0.5tpm[,c(1:15,19:21,28:30)]
head(FLB.control.count.data.0.5tpm.timepoints)
dim(FLB.control.count.data.0.5tpm.timepoints)


# so now I have a dataframe with counts for genes > 0.5 tpm in RNAi 
head(FLB.RNAi.count.data.0.5tpm)
# and for control only using the same timepoints as RNAi >0.5 tpm
head(FLB.control.count.data.0.5tpm.timepoints)


### DESeq2 ####

counts_for_DESeq2 <- merge(FLB.RNAi.count.data.0.5tpm, FLB.control.count.data.0.5tpm.timepoints, by = 0)
head(counts_for_DESeq2)
rownames(counts_for_DESeq2) <- counts_for_DESeq2[,1]
counts_for_DESeq2 <- counts_for_DESeq2[,-1]
head(counts_for_DESeq2)
dim(counts_for_DESeq2)

library(DESeq2)
colnames(counts_for_DESeq2)
timepoints <- c("FLB3_RNAi", "FLB7_RNAi", "FLB10_RNAi", "FLB13_RNAi", "FLB15_RNAi", "FLB19_RNAi","FLB26_RNAi",
                "FLB3_WT", "FLB7_WT", "FLB10_WT", "FLB13_WT", "FLB15_WT", "FLB19_WT","FLB26_WT")
CondVector <- rep(timepoints,each=3)
CondVector

sampleTable <- data.frame(row.names=colnames(counts), condition=as.factor(CondVector))
sampleTable

counts <- round(counts_for_DESeq2)
head(counts)

dds <- DESeqDataSetFromMatrix(countData = counts, colData=sampleTable, design=~condition)
dds

#make sure the FLB3_RNAi is used as the reference condition :
dds$condition <- relevel(dds$condition, "FLB3_RNAi")

dim(dds)

# saved a copy of dds just in case
dds_copy <- dds

# run DeSeq2
dds <- DESeq(dds)

# will need to have the bckCDS_1 set correctly to match old script which used this name for dds
#uses tutorial from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/ tutorial:
bckCDS_1 <- dds


# want to save the results of DESeq (i.e. bckCDS_1)
# all timepoints + tissues
timepoints

times_list <- c("3","7","10","13","15","19","26")


# get data each timepoint RNAi vs WT
for (i in times_list) { 
  
  bck_res <- results(bckCDS_1,contrast=c("condition",paste0("FLB",i,"_WT"),paste0("FLB",i,"_RNAi")))
  
  # sort results on padj
  ordered_res <- bck_res[order(bck_res$padj),]
  head(ordered_res)
  tail(ordered_res)
  ordered_res_na.rm <- na.omit(ordered_res)
  head(ordered_res_na.rm)
  tail(ordered_res_na.rm)
  
  # output ordered_res to csv
  write.csv(ordered_res_na.rm[ordered_res_na.rm$padj<0.05,],file=paste0(i,"DAA_WT_vs_RNAi_results.csv"))
  
  assign(paste("FLB",i,"DAA_WT_vs_RNAi_.na.rm",sep=""), ordered_res_na.rm[ordered_res_na.rm$padj<0.05,]) 
  
  print(paste("FLB",i,"DAA_WT_vs_RNAi_.na.rm",sep=""))
  dim(assign(paste("FLB",i,"DAA_WT_vs_RNAi_.na.rm",sep=""), ordered_res_na.rm[ordered_res_na.rm$padj<0.05,]) )
  
}


DE_genes <- data.frame(timepoint = numeric(), 
                       DE_genes_0.05 = numeric(), upreg_0.05 = numeric(), downreg_0.05 = numeric(),
                       DE_genes_0.01 = numeric(), upreg_0.01 = numeric(), downreg_0.01 = numeric(),
                       DE_genes_0.001 = numeric(), upreg_0.001 = numeric(), downreg_0.001 = numeric())
DE_genes

for (i in times_list) { 

head(get(paste("FLB",i,"DAA_WT_vs_RNAi_.na.rm",sep="")))
  my_data <- get(paste("FLB",i,"DAA_WT_vs_RNAi_.na.rm",sep=""))
  head(my_data)
  DE_genes_0.05 <- nrow(my_data)
  upreg_0.05 <- nrow(my_data[my_data$log2FoldChange > 1,])
  downreg_0.05 <- nrow(my_data[my_data$log2FoldChange < -1,])
  DE_genes_0.01 <- nrow(my_data[my_data$padj < 0.01,])
  upreg_0.01 <- nrow(my_data[my_data$padj < 0.01 & my_data$log2FoldChange > 1,])
  downreg_0.01 <- nrow(my_data[my_data$padj < 0.01 & my_data$log2FoldChange < -1,])
  DE_genes_0.001 <- nrow(my_data[my_data$padj < 0.001,])
  upreg_0.001 <- nrow(my_data[my_data$padj < 0.001 & my_data$log2FoldChange > 1,])
  downreg_0.001 <- nrow(my_data[my_data$padj < 0.001 & my_data$log2FoldChange < -1,])

  DE_genes <- rbind(DE_genes, list(timepoint = as.numeric(i), 
                                DE_genes_0.05 = DE_genes_0.05, upreg_0.05 = upreg_0.05, downreg_0.05 = downreg_0.05,
                                DE_genes_0.01 = DE_genes_0.01, upreg_0.01 = upreg_0.01, downreg_0.01 = downreg_0.01,
                                DE_genes_0.001 = DE_genes_0.001, upreg_0.001 = upreg_0.001, downreg_0.001 = downreg_0.001))
  
}

DE_genes

colnames(DE_genes) <- c("timepoint","DE_genes_0.05","upreg2fold_0.05", "downreg2fold_0.05",
                        "DE_genes_0.01","upreg2fold_0.01", "downreg2fold_0.01",
                        "DE_genes_0.001","upreg2fold_0.001", "downreg2fold_0.001")



### now want to do GO term enrichment for each timepoint for the three thresholds ####

#### read in information about lengths and GO terms #########
# read in GO terms
all_go <- read.csv("C:\\Users\\borrillp\\Documents\\Travel_May2018\\RefSeqv1.1_masked\\previous_work\\IWGSC_stress_GO.csv", sep=",")
head(all_go)
all_go <- all_go[,c(1,2)]
colnames(all_go) <- c("Gene", "GO_term")
head(all_go)
dim(all_go)

# convert from v1.0 to v1.1 # this is probably not the perfect thing to do so it might need to be re-done.....
head(gsub("01G", "02G", all_go$Gene))

all_go$Gene <- (gsub("01G", "02G", all_go$Gene))
head(all_go)
dim(all_go)

all_go_HC <- all_go[!grepl("LC", all_go$Gene),]
head(all_go_HC)
dim(all_go_HC)

length(unique(all_go_HC$Gene)) # number of HC genes with go terms before removing ones which don't match v1.0 to v1.1

# only keep genes which were >99 % ID > 90% coverage from v1.0 to v1.1 
genes_to_transfer <- read.csv(file="C:\\Users\\borrillp\\Documents\\Travel_May2018\\RefSeqv1.1_masked\\previous_work\\genes_to_transfer_qcov90_pident99_same_ID.csv")
head(genes_to_transfer)

all_go <- all_go[all_go$Gene %in% genes_to_transfer$gene_v1.1,]
head(all_go)
dim(all_go)

length(unique(all_go$Gene)) # number of genes with go terms

# select only genes which were used for DESeq2
head(counts_for_DESeq2)
dim(counts_for_DESeq2)
all_go <- subset(all_go, Gene %in% rownames(counts_for_DESeq2))
dim(all_go)

#create vector for gene_lengths

# need to get lengths of genes not of transcripts
lengths <- read.csv(file=paste0(data_dir,"control_timecourse_gene_lengths.csv"), header=T)
head(lengths)
colnames(lengths) <- c("gene", "length")
head(lengths)

t1 <- subset(lengths, gene %in% rownames(counts_for_DESeq2))
head(t1)
dim(t1)

# turn into a vector called gene.lens to use with GOSeq
gene.lens <- as.numeric(t1$length)
names(gene.lens) = t1$gene
head(gene.lens)
length(gene.lens)


####### Do GO term enrichment ####
out_dir <- "GO_enrichment\\"

assayed.genes <- as.vector(t1$gene)
length(assayed.genes)

library(goseq)
#i=1

GO_enriched <- data.frame(category = character(),over_represented_pvalue = numeric(),
                          under_represented_pvalue = numeric(), numDEInCat = numeric(), 
                          numInCat = numeric(), term = character(), ontology = character(),
                          over_rep_padj = numeric(), timepoint_threshold = character())
head(GO_enriched)


# do the GO term enrichment for each timepoint
for (i in times_list) { 
  
  head(get(paste("FLB",i,"DAA_WT_vs_RNAi_.na.rm",sep="")))
  my_data <- get(paste("FLB",i,"DAA_WT_vs_RNAi_.na.rm",sep=""))
  head(my_data)
  upreg2fold_0.01 <- (my_data[my_data$padj < 0.01 & my_data$log2FoldChange > 1,])
  downreg2fold_0.01 <- (my_data[my_data$padj < 0.01 & my_data$log2FoldChange < -1,])
  upreg2fold_0.001 <- (my_data[my_data$padj < 0.001 & my_data$log2FoldChange > 1,])
  downreg2fold_0.001 <- (my_data[my_data$padj < 0.001 & my_data$log2FoldChange < -1,])
  

for (j in c("upreg2fold_0.01","downreg2fold_0.01","upreg2fold_0.001","downreg2fold_0.001")){
  genes_for_GO <- (get(j))
 # head(genes_for_GO)

  #now do GO stats analysis on the genes expressed in each pattern compared to all genes expressed
  #create a named binary vector for genes where one means differentially expressed and 0 means not differentially expressed
  de.genes <- rownames(genes_for_GO)
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  head(gene.vector)
  #now carry out the GOseq analysis
  pwf = nullp(gene.vector, bias.data = gene.lens, plot.fit = TRUE)
  GO.wall = goseq(pwf, gene2cat = all_go)
  #this gave table with p-values...now correct for multiple testing using FDR
  enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]
  head(enriched.GO)
  # add new column with over represented GO terms padj
  GO.wall$over_rep_padj=p.adjust(GO.wall$over_represented_pvalue, method="BH")
  write.table(GO.wall[GO.wall$over_rep_padj <0.05,], file = paste0(out_dir, i,"DAA_WT_vs_RNAi_",j, "_GOseq.tsv", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)
 
  GO_enriched_timepoint <- GO.wall[GO.wall$over_rep_padj <0.05 & GO.wall$ontology == "BP",]
  head(GO_enriched_timepoint)
  
  # if no enriched GO terms don't add to dataframe
  if(nrow(GO_enriched_timepoint)>0) {
    
    GO_enriched_timepoint$timepoint_threshold <- paste0(i,"DAA_",j)
    GO_enriched <- rbind(GO_enriched,GO_enriched_timepoint)
  }
   
}
}


# now want to re-arrange the table to be in an easier to interpret format
head(GO_enriched)
library(dplyr)
# need to add new column saying GO1, GO2 etc for each pattern
#GO_enriched$rank=unlist(with(GO_enriched,tapply(over_rep_padj,pattern,rank)))

GO_enriched_ranked <- GO_enriched %>%
  group_by(timepoint_threshold) %>%
  mutate(subrank = rank(over_rep_padj,ties.method = "first"))

#check rank worked ok
head(data.frame(GO_enriched_ranked))
colnames(GO_enriched_ranked)

# now select columns I want to spread
GO_enriched_sel <- as.data.frame(GO_enriched_ranked[,c(6,9,10)])
head(GO_enriched_sel)

GO_data_spread <- spread(GO_enriched_sel, subrank, term)
GO_data_spread[1:5,1:5]

dim(GO_data_spread)
unique(GO_data_spread$timepoint_threshold)

# want to add number of genes in each DE group, and groups without BP GO terms enriched:
head(DE_genes_long)
tail(DE_genes_long)

# select only DE 2 fold padj 0.01 or 0.001 and 
# make a new column with the "timepoint_threshold" to facilitate merging
DE_genes_long_categories <- DE_genes_long %>%
  filter(threshold == "downreg2fold_0.001" |
           threshold == "upreg2fold_0.001" |
           threshold == "upreg2fold_0.01" |
           threshold == "downreg2fold_0.01" ) %>%
  mutate(timepoint_threshold = paste0(timepoint,"DAA_",threshold)) 

head(DE_genes_long_categories)

GO_data_spread[1:5,1:5]
merged_GO_data_spread <- merge(DE_genes_long_categories, GO_data_spread, by = "timepoint_threshold", all.x =T)
merged_GO_data_spread[1:10,1:5]

dim(merged_GO_data_spread)
unique(merged_GO_data_spread$timepoint_threshold)

#sort by timepoint
merged_GO_data_spread <- merged_GO_data_spread[order(merged_GO_data_spread$timepoint),]
merged_GO_data_spread[1:10,1:5]

# rename value column to num_DE_genes
colnames(merged_GO_data_spread)[1:10]
colnames(merged_GO_data_spread)[4] <- "num_DE_genes"

# separate into 0.01 and 0.001
merged_GO_data_spread_0.01 <- merged_GO_data_spread %>%
  filter(threshold == "upreg2fold_0.01" | 
           threshold == "downreg2fold_0.01")
merged_GO_data_spread_0.01[1:10,1:5]

merged_GO_data_spread_0.001 <- merged_GO_data_spread %>%
  filter(threshold == "upreg2fold_0.001" | 
           threshold == "downreg2fold_0.001")
merged_GO_data_spread_0.001[1:10,1:5]



write.csv(file=paste0(out_dir,"GO_enrichment_padj0.01.csv"), 
          merged_GO_data_spread_0.01, row.names = F)

write.csv(file=paste0(out_dir,"GO_enrichment_padj0.001.csv"), 
          merged_GO_data_spread_0.001, row.names = F)



