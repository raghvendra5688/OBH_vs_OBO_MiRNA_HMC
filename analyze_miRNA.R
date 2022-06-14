library(data.table)
library(ggplot2)
library(ggrepel)
library(readxl)
library(corrplot)
library(gplots)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(latex2exp)
library(stringr)
library(magick)

#Directory
setwd(".")
source("all_functions.R")

df_all_v1 <- read_excel("data/Fayaz_miRNA_FEB_2021_OBH_OBUH_NEW_20210214_153438_Export.xlsx",sheet = 1)
df_all_v1 <- as.data.frame(df_all_v1)
df_all_v1$Target.Name <- as.character(as.vector(df_all_v1$Target.Name))
df_all_v1$Biological.Group.Name <- as.character(as.vector(df_all_v1$Biological.Group.Name))
df_all_v1$Sample.Name <- as.character(as.vector(df_all_v1$Sample.Name))
df_all_v1$RQ <- as.numeric(as.vector(df_all_v1$RQ))
df_all_v1$`ΔCʀт Mean` <- as.numeric(as.vector(df_all_v1$`ΔCʀт Mean`))
df_all_v1$ΔΔCʀт <- as.numeric(as.vector(df_all_v1$ΔΔCʀт))
df_all_v1$P.Value <- as.numeric(as.vector(df_all_v1$P.Value))

#Get the subset comprising Obese unhealhty and Lean controls
subset_df_all_v1 <- df_all_v1[(df_all_v1$Biological.Group.Name=="OBUH" | df_all_v1$Biological.Group.Name=="OBH") & !is.na(df_all_v1$Biological.Group.Name),]

#Get all the patient information
patient_info <- read.table("data/common_patients.txt",header=TRUE,sep=",")
patient_info$Group <- as.character(as.vector(patient_info$Type))
patient_info$Sample_ID <- as.character(as.vector(patient_info$MiRNA_PatientID))

#Consider only patients belonging to the groups of interest
df_subset_v1 <- subset_df_all_v1[subset_df_all_v1$Sample.Name %in% patient_info$Sample_ID,]
unique_patients <- unique(df_subset_v1$Sample.Name)
for (patient in unique_patients)
{
  ids <- which(df_subset_v1$Sample.Name==patient)
  group <- patient_info[patient_info$Sample_ID==patient,]$Group
  df_subset_v1[ids,]$Biological.Group.Name=group
}
df_subset_v1 <- df_subset_v1[order(df_subset_v1$Sample.Name,df_subset_v1$Target.Name),]

#Get correct miRNA names
miRNA_info <- strsplit(df_subset_v1$Target.Name,"_")
miRNA_names <- NULL
for (i in 1:length(miRNA_info))
{
  temp <- miRNA_info[[i]][1]
  temp_name <- str_replace(string = temp, pattern = "hsa-", replacement = "")
  miRNA_names <- c(miRNA_names,temp_name)
}
df_subset_v1$miRNA_names <- miRNA_names

#Write output for obuh vs obh
out <- generate_plot_info(df_subset_v1, groupid = "OBUH", cutoff = 0.05)
g1 <- out[[1]]
ggsave(filename = "results/Volcano_plot_obuh_vs_obh_Fig1.pdf",plot = g1, device=pdf(), dpi=300, width=8, height=6, units="in")
dev.off()

#Write information about all microRNA
################################################################################
all_mirna_df <- out[[3]]
all_mirna_df$log2FC <- log2(all_mirna_df$RQ)
final_output_df <- all_mirna_df[, c(2,3,5,7)];
final_output_df <- final_output_df[order(final_output_df$P.Value),]
write.table(final_output_df,"results/all_microrna_info.csv",row.names=F, col.names=T,sep=",",quote=F)

#Get list of differentially expressed microRNAs
diff_miRNA_df <- all_mirna_df[all_mirna_df$P.Value<=0.05,]
diff_miRNA_df$ΔCʀт.Mean <- signif(diff_miRNA_df$ΔCʀт.Mean,digits=3)
diff_miRNA_df$RQ <- signif(diff_miRNA_df$RQ,digits=3)
diff_miRNA_df$P.Value <- signif(diff_miRNA_df$P.Value,digits=3)
write.table(diff_miRNA_df[,c(2:5)],"results/Diff_MicroRNAs.csv", row.names=F, col.names=T, sep=",",quote=F)

#Compare the Delta CRT of diff miRNA in OBUH and OBH samples
################################################################################
diff_df_subset_v1 <- df_subset_v1[df_subset_v1$miRNA_names %in% diff_miRNA_df$miRNA_names,]
unique_miRNAs <- unique(diff_df_subset_v1$miRNA_names)
miRNA_exp_matrix <- matrix(0, nrow=length(unique_miRNAs), ncol=2)
logpval_vector <- rep(0,length(unique_miRNAs))
logRQ_vector <- rep(0, length(unique_miRNAs))
for (i in 1:length(unique_miRNAs))
{
  miRNA_name <- unique_miRNAs[i]
  obuh_delta_crt_mean <- mean(diff_df_subset_v1[diff_df_subset_v1$Biological.Group.Name=="OBUH" & diff_df_subset_v1$miRNA_names==miRNA_name,]$`ΔCʀт Mean`)
  obuh_logpval_mean <- mean(-log10(diff_df_subset_v1[diff_df_subset_v1$Biological.Group.Name=="OBUH" & diff_df_subset_v1$miRNA_names==miRNA_name,]$P.Value))
  obuh_logRQ_mean <- mean(log2(diff_df_subset_v1[diff_df_subset_v1$Biological.Group.Name=="OBUH" & diff_df_subset_v1$miRNA_names==miRNA_name,]$RQ))  
  logpval_vector[i] <- obuh_logpval_mean
  logRQ_vector[i] <- obuh_logRQ_mean
  obh_delta_crt_mean <- mean(diff_df_subset_v1[diff_df_subset_v1$Biological.Group.Name=="OBH" & diff_df_subset_v1$miRNA_names==miRNA_name,]$`ΔCʀт Mean`)
  miRNA_exp_matrix[i,c(1,2)] <- c(obh_delta_crt_mean, obuh_delta_crt_mean)
}
colnames(miRNA_exp_matrix) <- c("OBO","OBM")
rownames(miRNA_exp_matrix) <- unique_miRNAs

groups <- c("OBO","OBM")
groups <- factor(groups, levels=c("OBO","OBM"))
col_fun1 <- colorRamp2(c(-2,0,2),c("blue","white","red"))
col_fun2 <- colorRamp2(c(1,2.5), c("blue","red"))

ht_fc <- make_heatmap_RQ(miRNA_exp_matrix, groups, col_fun1, col_fun2)
pdf(file="results/Diff_miRNA_mean_delta_crt_exp.pdf", height = 12, width=6, pointsize=11)
draw(ht_fc)
dev.off()

#Make the correlations of miRNA CRt expression with clinical traits
################################################################################
clinical_df <- read_excel("data/Clinical_data_common.xlsx",sheet = 1)
subset_clinical_df <- clinical_df[clinical_df$code %in% patient_info$PatientID,]
subset_clinical_df <- subset_clinical_df[order(subset_clinical_df$code),]
colnames(subset_clinical_df)[c(4,5,16)] <- c("Race","Gender","HBA1c")
patient_info <- patient_info[order(patient_info$PatientID),]
subset_clinical_df$patient_id <- patient_info$MiRNA_PatientID
features_of_interest <- c("patient_id","Race","Gender","Smoking","Age","HBA1c","HB","UREA","Creatinine","ALBUMIN","ALT","AST",
                          "CHOLESTROLE","TG","HDL","LDL","vitamin D","TSH","T4","vit B12")
rev_subset_clinical_df <- as.data.frame(subset_clinical_df[,features_of_interest])

#Build the ddCRt values matrix for miRNAs
unique_samples <- unique(diff_df_subset_v1$Sample.Name)
diff_df_subset <- as.data.frame(diff_df_subset_v1)
miRNA_crt_mean_exp_matrix <- matrix(NA, nrow=length(unique_samples), ncol=length(unique_miRNAs))
rownames(miRNA_crt_mean_exp_matrix) <- unique_samples
colnames(miRNA_crt_mean_exp_matrix) <- unique_miRNAs
for (i in 1:length(unique_samples))
{
  sample_name <- unique_samples[i]
  for (j in 1:length(unique_miRNAs))
  {
    miRNA_name <- unique_miRNAs[j]
    temp_df <- as.data.frame(diff_df_subset_v1[diff_df_subset_v1$Sample.Name==sample_name & grepl(miRNA_name,diff_df_subset$Target.Name),])
    temp_crt <- mean(as.numeric(as.vector(temp_df$Cʀт)),na.rm = T)
    miRNA_crt_mean_exp_matrix[sample_name, miRNA_name] <- temp_crt
  }
}

#remove diff miRNAs with no expression in over 10 patients
#################################################
no_of_NAs <- colSums(is.na(miRNA_crt_mean_exp_matrix))
remove_NAs_ids <- which(no_of_NAs>10)
rev_miRNA_crt_mean_exp_matrix <- miRNA_crt_mean_exp_matrix[,-remove_NAs_ids]
rev_miRNA_crt_mean_exp_matrix <- rev_miRNA_crt_mean_exp_matrix[order(rownames(rev_miRNA_crt_mean_exp_matrix)),]
rownames(rev_subset_clinical_df) <- rev_subset_clinical_df[,1]
rev_subset_clinical_df <- rev_subset_clinical_df[order(rownames(rev_subset_clinical_df)),]
final_subset_clinical_df <- rev_subset_clinical_df[rownames(rev_miRNA_crt_mean_exp_matrix),-c(1)]

#Build the correlation and pvalue matrix
correlation_matrix <- matrix(0, nrow=ncol(rev_miRNA_crt_mean_exp_matrix), ncol = ncol(final_subset_clinical_df))
pval_matrix <- matrix(0, nrow=ncol(rev_miRNA_crt_mean_exp_matrix), ncol = ncol(final_subset_clinical_df))
orig_miRNAs_of_interest <- colnames(rev_miRNA_crt_mean_exp_matrix)
clinical_targets_of_interest <- colnames(final_subset_clinical_df)
rownames(correlation_matrix) <- orig_miRNAs_of_interest
colnames(correlation_matrix) <- clinical_targets_of_interest
rownames(pval_matrix) <- orig_miRNAs_of_interest
colnames(pval_matrix) <- clinical_targets_of_interest
for (i in 1:ncol(rev_miRNA_crt_mean_exp_matrix))
{
  miRNA_name <- orig_miRNAs_of_interest[i]
  for (j in 1:ncol(final_subset_clinical_df))
  {
    clinical_target_name <- clinical_targets_of_interest[j]
    cor_info <- cor.test(x=-rev_miRNA_crt_mean_exp_matrix[,miRNA_name],y=final_subset_clinical_df[,clinical_target_name],method = "pearson",
                       na.action=na.omit)
    cor_p_value <- cor_info$p.value
    cor_value <- as.numeric(cor_info$estimate)
    correlation_matrix[miRNA_name, clinical_target_name] <- cor_value
    pval_matrix[miRNA_name, clinical_target_name] <- cor_p_value
  }
}

metabolic_features <- c("HBA1c","Creatinine","CHOLESTROLE","LDL","HDL")
rev_correlation_matrix <- correlation_matrix[,metabolic_features]
rev_pval_matrix <- pval_matrix[,metabolic_features]
miRNAs_of_interest <- names(which(rowSums(rev_pval_matrix<0.05)>0))
final_correlation_matrix <- rev_correlation_matrix[miRNAs_of_interest,]
final_pval_matrix <- rev_pval_matrix[miRNAs_of_interest,]

pdf("results/miRNA_clinical_trait_correlation_ddCRT.pdf",height=6, width =4)
corrplot(corr=final_correlation_matrix,
         p.mat = final_pval_matrix,
         type="full", insig="pch", sig.level =.05, pch.cex = 1.0, col=bluered(100))
dev.off()

#Keep only 3 digits for significance
for(i in 1:nrow(final_correlation_matrix))
{
  for(j in 1:ncol(final_correlation_matrix))
  {
    final_correlation_matrix[i,j] <- signif(final_correlation_matrix[i,j],3)
    final_pval_matrix[i,j] <- signif(final_pval_matrix[i,j],3)
  }
}
write.table(final_correlation_matrix, file="results/miRNA_clinical_trait_correlation_matrix_ddCRT.csv",
            row.names=T, col.names=T, quote=F, sep="\t")
write.table(final_pval_matrix, file="results/miRNA_clinical_trait_pval_matrix_ddCRT.csv",
            row.names=T, col.names=T, quote=F, sep="\t")

#Make scatter plot for HBA1c etc
rownames(patient_info) <- patient_info$MiRNA_PatientID
subset_patient_info <- patient_info[rownames(rev_miRNA_crt_mean_exp_matrix),]
revised_subset_clinical_df <- rev_subset_clinical_df[rownames(rev_miRNA_crt_mean_exp_matrix),]
g_hba1c <- make_scatter_plot(subset_patient_info, rev_miRNA_crt_mean_exp_matrix, revised_subset_clinical_df, pval_matrix, trait="HBA1c")
ggsave(filename=paste0("results/HBA1c_significant_correlated_diff_miRNAs_ddCRT.pdf"), plot = g_hba1c, device = pdf(), height=5 ,width = 10,units="in")
dev.off()

g_creatinine <- make_scatter_plot(subset_patient_info, rev_miRNA_crt_mean_exp_matrix, revised_subset_clinical_df, pval_matrix, trait="Creatinine")
ggsave(filename=paste0("results/Creatinine_significant_correlated_diff_miRNAs_ddCRT.pdf"),plot = g_creatinine, device = pdf(), height=5,width = 10,units="in")
dev.off()

g_cholesterol <- make_scatter_plot(subset_patient_info, rev_miRNA_crt_mean_exp_matrix, revised_subset_clinical_df, pval_matrix, trait="CHOLESTROLE")
ggsave(filename=paste0("results/Cholesterol_significant_correlated_diff_miRNAs_ddCRT.pdf"),plot = g_cholesterol, device = pdf(), height=5,width = 10,units="in")
dev.off()

g_hdl <- make_scatter_plot(subset_patient_info, rev_miRNA_crt_mean_exp_matrix, revised_subset_clinical_df, pval_matrix, trait="HDL")
ggsave(filename=paste0("results/HDL_significant_correlated_diff_miRNAs_ddCRT.pdf"),plot = g_hdl, device = pdf(), height=5,width = 3.5,units="in")
dev.off()

g_ldl <- make_scatter_plot(subset_patient_info, rev_miRNA_crt_mean_exp_matrix, revised_subset_clinical_df, pval_matrix, trait="LDL")
ggsave(filename=paste0("results/LDL_significant_correlated_diff_miRNAs_ddCRT.pdf"),plot = g_ldl, device = pdf(), height=5,width = 3.5,units="in")
dev.off()

#Keep only miRNAs with significant correlations
miRNA_with_significant_correlations <- union(union(names(which(pval_matrix[,"HBA1c"]<0.05)),
                                             names(which(pval_matrix[,"Creatinine"]<0.05))),
                                             names(which(pval_matrix[,"CHOLESTROLE"]<0.05)))

#Get the diff miRNA-mRNA gene with strong evidence for interaction from MirDIP database
#################################################################################
all_miRNA_interactions_df <- read_excel("data/mirDIP_diff_miRNA_mRNA_interactions.xlsx",sheet=1)

#Consider only those interactions between diff miRNAs and mRNAs with >=10 sources 
all_miRNA_interactions_high_score_df <- all_miRNA_interactions_df[all_miRNA_interactions_df$`Number of Sources`>=10 & all_miRNA_interactions_df$`Integrated Score`>0.75,]
edgelist_df <- all_miRNA_interactions_high_score_df[,c("MicroRNA","Gene Symbol","Integrated Score")]
colnames(edgelist_df) <- c("Source","Target","Weight")
edgelist_df$Source <- str_replace(edgelist_df$Source, pattern="hsa-", replacement = "")
subset_edgelist_df <- edgelist_df[edgelist_df$Source %in% miRNA_with_significant_correlations,]
write.table(subset_edgelist_df, file="results/Diff_miRNA_mRNA_interactions_to_plot_ddCRT.csv",row.names=F, col.names=T, sep=",", quote=F)

all_nodes <- union(unique(subset_edgelist_df$Source),unique(subset_edgelist_df$Target))
node_type <- rep(0, length(all_nodes))
node_type[which(all_nodes %in% diff_miRNA_df[diff_miRNA_df$log2FC>0,]$miRNA_names)] <- 1
node_type[which(all_nodes %in% diff_miRNA_df[diff_miRNA_df$log2FC<0,]$miRNA_names)] <- 2
nodetable_df <- data.frame(Nodes = all_nodes, Type = node_type )
write.table(nodetable_df, file="results/All_miRNA_mRNA_nodetable_ddCRT.csv", row.names=F, col.names=T, sep=",", quote=F)

rev_full_edgelist <- all_miRNA_interactions_high_score_df[all_miRNA_interactions_high_score_df$MicroRNA %in% miRNA_with_significant_correlations &
                                                            all_miRNA_interactions_high_score_df$`Number of Sources`>=10 &
                                                            all_miRNA_interactions_high_score_df$`Integrated Score`>0.75,]
write.table(rev_full_edgelist, file="results/Diff_miRNA_mRNA_interactions_full_list_ddCRT.csv",row.names=F, col.names=T, sep="\t", quote=F)

#See the impact of sex and smoking 
################################################################################
combined_df <- as.data.frame(cbind(subset_patient_info$Group, final_subset_clinical_df[,c(1:4)], -rev_miRNA_crt_mean_exp_matrix))
colnames(combined_df)[1] <- "Group"
combined_df$Group <- as.factor(as.vector(combined_df$Group))
combined_df$Race <- as.factor(as.vector(combined_df$Race))
combined_df$Gender <- as.factor(as.vector(combined_df$Gender))
combined_df$Smoking <- as.factor(as.vector(combined_df$Smoking))

covariates <- c(3,4)
delta_crt_matrix <- matrix(0, nrow=length(miRNAs_of_interest), ncol=2)
pval_matrix <- matrix(1, nrow=length(miRNAs_of_interest), ncol=2)
rownames(delta_crt_matrix) <- miRNAs_of_interest
rownames(pval_matrix) <- miRNAs_of_interest
for (i in 1:length(miRNAs_of_interest))
{
  for (j in 1:length(covariates))
  {
    temp_wt <- t.test(x=combined_df[combined_df[,covariates[j]]==0,miRNAs_of_interest[i]],y=combined_df[combined_df[,covariates[j]]==1,miRNAs_of_interest[i]], exact = F)
    delta_crt_matrix[miRNAs_of_interest[i],j] <- as.numeric(temp_wt$estimate[2]-temp_wt$estimate[1])
    pval_matrix[miRNAs_of_interest[i],j] <- temp_wt$p.value
  }
}
colnames(pval_matrix) <- c("Gender","Smoking")
colnames(delta_crt_matrix) <- c("Gender","Smoking")

pdf("results/miRNA_gender_smoking_trait_ttest_ddCRT.pdf",height=6, width =4)
ht_dcrt <- make_heatmap_dcrt(delta_crt_matrix, pval_matrix)
print(ht_dcrt)
dev.off()

g_gender <- make_boxplot_plot(patient_info, rev_miRNA_crt_mean_exp_matrix, rev_subset_clinical_df, pval_matrix, trait="Gender")
ggsave(filename="results/Gender_significant_ttest_diff_miRNAs_ddCRT.pdf", plot = g_gender, device = pdf(), height=5,width = 7,units="in" )
dev.off()
g_smoking <- make_boxplot_plot(patient_info, rev_miRNA_crt_mean_exp_matrix, rev_subset_clinical_df, pval_matrix, trait="Smoking")
ggsave(filename="results/Smoking_significant_ttest_diff_miRNAs_ddCRT.pdf", plot = g_smoking, device = pdf(), height=5,width = 3.5,units="in" )
dev.off()
