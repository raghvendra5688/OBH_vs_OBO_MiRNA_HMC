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

setwd(".")

#Make the function to get volcano plot and differentially expressed info
generate_plot_info <- function(df_all,groupid,cutoff)
{
  #Get OBUH vs OBH fold-change data using RQ
  a_vs_b_full <- df_all[df_all$Biological.Group.Name==groupid,]
  min_RQ_log2 <- round(log2(min(a_vs_b_full$RQ)),2)
  max_RQ_log2 <- round(log2(max(a_vs_b_full$RQ)),2)
  max_Pvalue <- round(-log10(min(a_vs_b_full$P.Value)),2)
  min_Pvalue <- round(-log10(max(a_vs_b_full$P.Value)),2)
  unique_miRNAs <- unique(a_vs_b_full$miRNA_names)
  
  a_vs_b <- NULL
  for (i in 1:length(unique_miRNAs))
  {
    miRNA <- unique_miRNAs[i]
    temp_df <- a_vs_b_full[a_vs_b_full$miRNA_names==miRNA,]
    temp <- cbind(unique(temp_df$Biological.Group.Name),unique(temp_df$miRNA_names),
                  mean(temp_df$ΔΔCʀт),mean(temp_df$RQ),mean(temp_df$P.Value))
    a_vs_b <- rbind(a_vs_b,temp)
  }
  a_vs_b <- as.data.frame(a_vs_b)
  colnames(a_vs_b) <- c("Biological_Group","miRNA_names","ΔCʀт.Mean","RQ","P.Value")
  a_vs_b$RQ <- as.numeric(as.vector(a_vs_b$RQ))
  a_vs_b$P.Value <- as.numeric(as.vector(a_vs_b$P.Value))
  a_vs_b$ΔCʀт.Mean <- as.numeric(as.vector(a_vs_b$ΔCʀт.Mean))
  colors <- rep("black",nrow(a_vs_b))
  colors[log2(a_vs_b$RQ)>=2 & a_vs_b$P.Value < 0.05] <- "red"
  colors[log2(a_vs_b$RQ)<=-2 & a_vs_b$P.Value < 0.05] <- "blue"
  a_vs_b$colors <- colors
  
  #Filtered data for text annotation
  annotation_df <- a_vs_b[-log10(a_vs_b$P.Value)>=-log10(cutoff) & abs(log2(a_vs_b$RQ))>=2.0,]
  annotation_df$log2FC <- log2(annotation_df$RQ)
  
  #Make plot for log2-fold change vs P.Value
  g <- ggplot(data=a_vs_b, aes(x=log2(RQ), y=-log10(P.Value))) +
    geom_point(size=4 ,aes(color=colors)) +
    theme(legend.position = "none") +
    ggtitle(paste0("Volcano Plot: log2(RQ) vs P-Value for OBM vs OBO")) +
    geom_hline(yintercept = -log10(cutoff), col="blue") +
    geom_vline(xintercept = 2, col="black") +
    geom_vline(xintercept = -2, col="black") +
    xlab("log2(RQ)") + ylab("-log10(P-Value)")  +
    scale_color_manual(values=c("black","blue","red")) +
    ylim(c(0,4.0))+
    xlim(c(-6.0,6.0))+
    #geom_label_repel(data=annotation_df, aes(label = miRNA_names, x=log2(RQ), y=-log10(P.Value)),
    #                 box.padding   = 0.25, max.overlaps = 50,
    #                 point.padding = 0.25,
    #                 segment.color = 'grey50') +
    theme_light() +
    theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "bold"),
          axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "bold"),
          axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "bold"),
          axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "bold"),
          legend.position = "none",
          text = element_text(size=20, face="bold"))
  return(list(g,annotation_df,a_vs_b))
}


make_heatmap_RQ <- function(miRNA_exp_matrix, groups, col_fun1, col_fun2)
{
  ht_fc = Heatmap(-miRNA_exp_matrix, col = colorRamp2(c(-5,0,5),c("blue","white","red")),
                  heatmap_legend_param = list(title = TeX("-$\\Delta CR_{T}$ Mean")),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height,
                              gp = gpar(col = "grey", fill = NA))
                    
                  },
                  width = unit(4, "cm"),
                  height = unit(24, "cm"),
                  cluster_columns = F,
                  cluster_rows = T,
                  row_names_centered = F,
                  show_row_names = T,
                  row_labels = rownames(miRNA_exp_matrix),
                  row_names_max_width = max_text_width(rownames(miRNA_exp_matrix), gp = gpar(fontsize = 11)),
                  right_annotation = rowAnnotation(logRQ = logRQ_vector,
                                                   logPval = logpval_vector,
                                                   simple_anno_size = unit(4,"mm"),
                                                   annotation_name_gp = gpar(fontsize=11, family="sans", col="black"),
                                                   col = list(logRQ = col_fun1,
                                                              logPval = col_fun2),
                                                   annotation_legend_param = list(direction="vertical"),
                                                   show_legend=T),                
                  column_split = groups,
                  show_column_names = F,
                  raster_quality = 2,
                  row_names_side = "right",
                  use_raster = T,
                  column_gap = unit(3, "mm"),
                  border_gp = gpar(col = "black", lty = 1),
                  border = T
  )
  return(ht_fc)
}

#Make individual correlation scatter plot for clinical traits at significant correlation with at least 3 different miRNAs: HBA1c, HB, Creatinine, Cholestrole
make_scatter_plot <- function(patient_info, rev_miRNA_crt_mean_exp_matrix, rev_subset_clinical_df, pval_matrix, trait="HBA1c")
{
  rev_patient_info <- patient_info[order(patient_info$MiRNA_PatientID),]
  miRNA_names <- names(which(pval_matrix[,trait]<0.05))
  temp_clinical_df <- as.data.frame(cbind(rev_subset_clinical_df[,c("patient_id",trait)],rev_patient_info$Group,-rev_miRNA_crt_mean_exp_matrix[,names(which(pval_matrix[,trait]<0.05))]))
  colnames(temp_clinical_df) <- c("patient_id",trait,"Group",miRNA_names)
  temp_clinical_df <- reshape2::melt(temp_clinical_df, id=c("patient_id",trait,"Group"))
  rev_temp_clinical_df <- temp_clinical_df[complete.cases(temp_clinical_df),]
  variables <- unique(rev_temp_clinical_df$variable)
  g <- ggplot(data = rev_temp_clinical_df, aes_string(x=trait, y="value")) + facet_wrap(~variable, nrow = 2, ncol=3) + geom_point(aes(color=Group)) +
    stat_cor(method="pearson", p.accuracy = 0.001, r.accuracy = 0.01)+ylab(TeX("-$CR_{T}$ value"))+
    geom_smooth(method="lm",col="black") + scale_color_discrete(labels=c("OBO","OBM"),type=c("blue","red")) +
    theme_bw() + theme(text = element_text(size=16),
                       axis.text.x = element_text(angle=0, hjust=1, size=16),
                       axis.text.y = element_text(size=16))
  
  return(g)
}

make_heatmap_dcrt <- function(delta_crt_matrix, pval_matrix)
{
  ht_dcrt = Heatmap(delta_crt_matrix, col = colorRamp2(c(-1,0,1),c("blue","white","red")),
                  heatmap_legend_param = list(title = TeX("$\\Delta -CR_{T}$")),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height,
                              gp = gpar(col = "grey", fill = NA))
                    if (pval_matrix[i,j]<=0.05)
                    {
                      grid.text(sprintf("*"), x, y, gp = gpar(fontsize = 11))
                    }
                    
                  },
                  width = unit(4, "cm"),
                  height = unit(10, "cm"),
                  cluster_columns = F,
                  cluster_rows = F,
                  row_names_centered = F,
                  show_row_names = T,
                  row_labels = rownames(delta_crt_matrix),
                  row_names_max_width = max_text_width(rownames(delta_crt_matrix), gp = gpar(fontsize = 11)),
                  column_split = c("Gender","Smoking"),
                  show_column_names = F,
                  raster_quality = 2,
                  row_names_side = "right",
                  use_raster = T,
                  column_gap = unit(3, "mm"),
                  border_gp = gpar(col = "black", lty = 1),
                  border = T
  )
  return(ht_dcrt)
}

make_boxplot_plot <- function(patient_info, rev_miRNA_crt_mean_exp_matrix, rev_subset_clinical_df, pval_matrix, trait="Gender")
{
  rev_patient_info <- patient_info[order(patient_info$MiRNA_PatientID),]
  miRNA_names <- names(which(pval_matrix[,trait]<0.05))
  temp_clinical_df <- as.data.frame(cbind(rev_subset_clinical_df[,c("patient_id",trait)],rev_patient_info$Group,-rev_miRNA_crt_mean_exp_matrix[,names(which(pval_matrix[,trait]<0.05))]))
  colnames(temp_clinical_df) <- c("patient_id",trait,"Group",miRNA_names)
  temp_clinical_df <- reshape2::melt(temp_clinical_df, id=c("patient_id",trait,"Group"))
  rev_temp_clinical_df <- temp_clinical_df[complete.cases(temp_clinical_df),]
  rev_temp_clinical_df[,trait] <- as.factor(as.vector(rev_temp_clinical_df[,trait]))
  variables <- unique(rev_temp_clinical_df$variable)
  g <- ggplot(data = rev_temp_clinical_df, aes_string(x=trait, y="value")) + facet_wrap(~variable, nrow = 2, ncol=3) + 
    #stat_cor(method="pearson", p.accuracy = 0.001, r.accuracy = 0.01)+
    geom_boxplot()+stat_compare_means(method = "t.test")+geom_jitter(aes(color=Group), width=0.25) +
    ylab(TeX("-$CR_{T}$ value"))+
    #geom_smooth(method="lm",col="black") + 
    scale_color_discrete(labels=c("OBO","OBM"),type=c("blue","red")) +
    theme_bw() + theme(text = element_text(size=16),
                       axis.text.x = element_text(angle=0, hjust=1, size=16),
                       axis.text.y = element_text(size=16))
  
  return(g)
}