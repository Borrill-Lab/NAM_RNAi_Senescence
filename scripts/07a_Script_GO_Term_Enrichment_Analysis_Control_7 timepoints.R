# Aim is to look at DE genes found by ImpulseDE2 and gradient tool output, and group genes into clusters which have the same expression pattern and do GO enrichment on clusters
# 20.7.2021
# Tayyaba Andleeb


library(dplyr)
library(tidyr)
library(ggplot2)

out_dir <- "C:/Users/TXA013/3D Objects/Borril Lab/Go enrichment term analysis/Go_Term_Enrichment_Analysis_Control_7 timepoints/"


# load genes DE in ImpulseDE2 + gradient tool Control_7timepoints
Control_DE <- read.csv(file="C:/Users/TXA013/3D Objects/Borril Lab/ImpulseDE2/ImpulseDE2_Control_7timepoints/merged_impuseDE2_padj0.001_gradient_tool_tpm_normalised_7timepoints.csv")
head(Control_DE)
dim(Control_DE)

# now keep only DE in gradient tool and impulseDE2
Control_DE <- Control_DE[!is.na(Control_DE$pattern) & !is.na(Control_DE$padj),]
dim(Control_DE)
head(Control_DE)

# how many genes are in each pattern (not much point in keeping patterns with very few genes - should group them together- but how to decide that?)
library(dplyr)
library(magrittr)
number_per_pattern <- Control_DE %>%
  group_by(pattern) %>%
  summarize(n= length(pattern))

write.csv(file=paste0(out_dir,"Control_7 timepoints_number_per_pattern.csv"), number_per_pattern, row.names = F)

head(number_per_pattern)
nrow(number_per_pattern) # number of patterns
nrow(as.data.frame(number_per_pattern)[number_per_pattern$n > 10,]) # number of patterns with >10 genes
nrow(as.data.frame(number_per_pattern)[number_per_pattern$n > 20,]) # number of patterns with >10 genes


### MANUALLY curated patterns in excel ###
# load grouped patterns
grouped_patterns <- read.csv(file=paste0(out_dir,"Control_7 timepoints_number_per_pattern_grouped.csv"),)
head(grouped_patterns)

# add in grouped pattern info to RNAi_DE
Control_DE <- merge(Control_DE, grouped_patterns, by.x = "pattern", by.y="pattern", all.x=T)
head(Control_DE)
dim(Control_DE)
sum(complete.cases(Control_DE))# check all rows have data

# save WT_7time points_DE with grouped patterns
write.csv(file=paste0(out_dir,"Control_7 timepoints_DE_with_grouped_pattern.csv"), Control_DE, row.names = F)



### do GO enrichment on clusters

#### read in information about lengths and GO terms #########
# read in GO terms

all_go <- read.csv("C:/Users/TXA013/3D Objects/Borril Lab/Go enrichment term analysis/Copy of IWGSC_stress_GO.csv") 
head(all_go)
all_go <- all_go[,c(1,2)]
colnames(all_go) <- c("Gene", "GO_term")
head(all_go)
dim(all_go)

# convert from v1.0 to v1.1 
head(gsub("01G", "02G", all_go$Gene))

all_go$Gene <- (gsub("01G", "02G", all_go$Gene))
head(all_go)
dim(all_go)

all_go_HC <- all_go[!grepl("LC", all_go$Gene),]
head(all_go_HC)
dim(all_go_HC)

length(unique(all_go_HC$Gene)) # number of HC genes with go terms before removing ones which don't match v1.0 to v1.1


# only keep genes which were >99 % ID > 90% coverage from v1.0 to v1.1 # avoids erroneous transfer of GO terms from v1.0 to v1.1
genes_to_transfer <- read.csv(file="C:/Users/TXA013/3D Objects/Borril Lab/Go enrichment term analysis/Copy of genes_to_transfer_qcov90_pident99_same_ID.csv")
head(genes_to_transfer)

all_go <- all_go[all_go$Gene %in% genes_to_transfer$gene_v1.1,]
head(all_go)
dim(all_go)

length(unique(all_go$Gene)) # number of genes with go terms

# select only genes which were used for gradient tool (i.e. remove non-expressed genes) 
expr_genes <- read.csv("C:/Users/TXA013/3D Objects/Borril Lab/Go enrichment term analysis/Go_Term_Enrichment_Analysis_Control_7 timepoints/control_7timepoints_tpm_FLB_0.5tpm_input_gradient_tool.csv")
head(expr_genes)
all_go <- subset(all_go, Gene %in% expr_genes$X)
dim(all_go)
length(unique(expr_genes$X)) # number of genes expressed

length(unique(all_go$Gene)) #  number of genes with go terms which were expressed

#create vector for gene_lengths

# need to get lengths of genes not of transcripts
lengths <- read.csv(file="C:/Users/TXA013/3D Objects/Borril Lab/Go enrichment term analysis/Copy of control_timecourse_gene_lengths.csv", header=T)
head(lengths)
colnames(lengths) <- c("gene", "length")
head(lengths)

t1 <- subset(lengths, gene %in% expr_genes$X)
head(t1)
dim(t1)

# turn into a vector called gene.lens to use with GOSeq
gene.lens <- as.numeric(t1$length)
names(gene.lens) = t1$gene
head(gene.lens)
length(gene.lens)

####### Do GO term enrichment and plot graph for each group ####

assayed.genes <- as.vector(t1$gene)
length(assayed.genes)
head(assayed.genes)
head(Control_DE)
library(BiasedUrn)
library(geneLenDataBase)
library(goseq)
library(dplyr)
library(tidyr)
library(org.Mm.eg.db)
#i=1

# do the GO term enrichment for each cluster

GO_enriched <- data.frame(category = character(),over_represented_pvalue = numeric(),
                          under_represented_pvalue = numeric(), numDEInCat = numeric(), 
                          numInCat = numeric(), term = character(), ontology = character(),
                          over_rep_padj = numeric(),grouped_pattern=character())
head(GO_enriched)

for (i in (unique(Control_DE$grouped_pattern))){
  
  #now do GO stats analysis on the genes expressed in each pattern compared to all genes expressed
  #create a named binary vector for genes where one means differentially expressed and 0 means not differentially expressed
  de.genes <- (Control_DE[Control_DE$grouped_pattern == i,"X"])
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  head(gene.vector)
  #now carry out the GOseq analysis
  pwf = nullp(gene.vector, bias.data = gene.lens, plot.fit = TRUE)
  GO.wall = goseq(pwf, gene2cat = all_go)
  # add new column with over represented GO terms padj
  GO.wall$over_rep_padj=p.adjust(GO.wall$over_represented_pvalue, method="BH")
  write.table(GO.wall[GO.wall$over_rep_padj <0.05,], file = paste0(out_dir, "cluster_", i, "_GOseq.tsv", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)
  
  GO_enriched_cluster <- GO.wall[GO.wall$over_rep_padj <0.05 & GO.wall$ontology == "BP",]
  head(GO_enriched_cluster)
  
  # if no enriched GO terms don't add to dataframe
  if(nrow(GO_enriched_cluster)>0) {
    
    GO_enriched_cluster$pattern <- i
    GO_enriched <- rbind(GO_enriched,GO_enriched_cluster)
  }
  #GO_enriched
  
}

# now want to re-arrange the table to be in an easier to interpret format
head(GO_enriched)
# need to add new column saying GO1, GO2 etc for each pattern
#GO_enriched$rank=unlist(with(GO_enriched,tapply(over_rep_padj,pattern,rank)))

GO_enriched_ranked <- GO_enriched %>%
  group_by(pattern) %>%
  mutate(subrank = rank(over_rep_padj,ties.method = "first"))

#head(data.frame(GO_enriched_ranked),200)
colnames(GO_enriched_ranked)

# now select columns I want to spread
GO_enriched_sel <- as.data.frame(GO_enriched_ranked[,c(6,9,10)])
head(GO_enriched_sel)

GO_data_spread <- spread(GO_enriched_sel, subrank, term)
GO_data_spread[1:5,1:5]



# want to add number of genes in each pattern, and patterns without BP GO terms enriched:dim(GO_data_spread)
head(grouped_patterns)
grouped_patterns %>%
  filter(grouped_pattern == "D03")

df_number_per_pattern <- grouped_patterns %>%
  group_by(grouped_pattern) %>%
  summarise(num_genes = sum(n))
head(df_number_per_pattern)

GO_data_spread[1:5,1:5]
merged_GO_data_spread <- merge(df_number_per_pattern, GO_data_spread, by.x = "grouped_pattern",by.y="pattern", all.x =T)
merged_GO_data_spread[1:10,1:5]

dim(merged_GO_data_spread)

# want to add column with the most common pattern of expression
head(grouped_patterns)
representative_pattern <- grouped_patterns %>% # pick the most common pattern for each grouped_pattern
  group_by(grouped_pattern) %>%
  summarise(num_genes_in_representative_pattern = max(n), 
            representative_pattern = pattern[which.max(n)])


merged_GO_data_spread_final <- merge(representative_pattern,merged_GO_data_spread, by.x = "grouped_pattern",by.y="grouped_pattern", all.x =T)
merged_GO_data_spread_final[1:10,1:20]

ncol(merged_GO_data_spread_final)
merged_GO_data_spread_final[1:5,139:145]

merged_GO_data_spread_final <- merged_GO_data_spread_final %>%
  dplyr::select(grouped_pattern,num_genes,representative_pattern,num_genes_in_representative_pattern,`1`:`139`)

# want to add columns with the chnage at each timepoint
merged_GO_data_spread_final <- merged_GO_data_spread_final %>%
  separate(representative_pattern, c("3D", "7D", "10D", "13D", "15D", "19D","26D"), "_", remove=F)
merged_GO_data_spread_final[1:7,1:20]
write.csv(file=paste0(out_dir,"GO_enrichment_per_cluster_control_7 timepoints.csv"), merged_GO_data_spread_final, row.names = F)