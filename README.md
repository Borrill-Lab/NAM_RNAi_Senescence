# NAM_RNAi_Senescence

This repository contains scripts and data used in the manuscript:

“Wheat NAM genes regulate the majority of early monocarpic senescence transcriptional changes including nitrogen remobilisation genes” by Tayyaba Andleeb and Philippa Borrill. 

# Data
Data files used are in <a href="https://github.com/Borrill-Lab/NAM_RNAi_Senescence/tree/main/data_files">data_files</a>.

v1.0 GO terms can be found here: <a href="https://github.com/Borrill-Lab/WheatFlagLeafSenescence/blob/master/data/IWGSC_stress_GO.csv"> IWGSC_stress_GO</a>

List of transcripts for which GO terms were transferred from the v1.0 to v1.1 gene annotation can be found here: <a href="https://github.com/Borrill-Lab/WheatFlagLeafSenescence/blob/master/data/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt"> transcripts v1.0 to v1.1 </a>

# Scripts

All scripts used are in <a href="https://github.com/Borrill-Lab/NAM_RNAi_Senescence/tree/main/scripts">scripts</a>.

<i>Mapping scripts: </i>

01a_kallisto_control.pl
Mapping of control (WT) samples to RefSeqv1.1 transcriptome.

01b_kallisto_RNAi.pl
Mapping of NAM RNAi samples to RefSeqv1.1 transcriptome.

02_tximport_summarise_counts_tpm_per_gene
Merge tpm and counts into a single table for all samples

03_WT_vs_RNAi_DESeq2_goseq
Identify genes differentially expressed between WT and RNAi at each timepoint using DESeq2. Carry out GO term enrichment analysis with goseq.

04a_Script_prep_data_control_7_timepoint_gradient_tool_cyverse.R
Prepare data for gradient tool on cyverse to look for genes differentially expressed across time in RNA-seq timecourse at seven timepoints in WT.

04b_Script_prep_data_RNAi_gradient_tool_cyverse.R
Prepare data for gradient tool on cyverse to look for genes differentially expressed across time in RNA-seq timecourse in RNAi.

05a_Script_gradient_tool_cluster_genes_control_7timespoints.R
Analyse gradient tool output and group genes into clusters which have the same expression pattern changes across seven timepoints in WT.

05b_Script_gradient tool_cluster_genes_RNAi.R
Analyse gradient tool output and group genes into clusters which have the same expression pattern changes across seven timepoints in RNAi.

06a_Script_ImpulseDE2_control_timecourse_7_timepoints.R
Run ImpulseDE2 on control timecourse for 7 timepoints and compare ImpulseDE2 results to gradient tool results in WT.

06b_Script_ImpulseDE2_RNAi.R
Run ImpulseDE2 on RNAi for 7 timepoints and compare ImpulseDE2 results to gradient tool results in RNAi.

07a_Script_GO_Term_Enrichment_Analysis_Control_7 timepoints.R
Identify DE genes found by ImpulseDE2 and gradient tool, and group genes into clusters which have the same expression pattern and do GO enrichment on clusters in WT. 

07b_Script_GO_Term_Enrichment_Analysis_RNAi.R
Identify DE genes found by ImpulseDE2 and gradient tool, and group genes into clusters which have the same expression pattern and do GO enrichment on clusters in RNAi.

08a_Script_Prepare_FLB_Avg_WT_7_timepoints.R
Calculate average TPM value and MaxTPM to select N responsive genes for 7 timepoints in WT.

08b_Script_Prepare_FLB_Avg_RNAi.R
Calculate average TPM value and MaxTPM to select N responsive genes for 7 timepoints in RNAi.
