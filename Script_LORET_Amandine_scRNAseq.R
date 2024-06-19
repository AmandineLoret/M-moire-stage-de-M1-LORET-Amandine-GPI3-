#Chargement des librairies
library(Seurat)
library(ggplot2)
library(dplyr)
library(plyr)
library(dittoSeq)
library(ggrepel)
library(enrichR)
library(WGCNA)
library(hdWGCNA)
library(patchwork)
library(cowplot)
library(tidyverse)
library(tidyr)
library(devtools)
library(presto)
library(ggpubr)

#Chargement de l'objet seurat clusterisé
pbmc_integre <- readRDS("/LAB-DATA/BiRD/shares/CRTI/EQ5/SC_pbmc/210928_pbmcEXCintegrated.rds")

#Création d'une colonne pour les conditions
pbmc_integre$condition <- NA
intersection <- data.frame(HTO = c("HTO1-TotalSeqB", "HTO2-TotalSeqB", "HTO3-TotalSeqB"), Condition = c("MS", "Ag_MS", "HV"))
for (i in 1:nrow(intersection)) {
  hto <- intersection$HTO[i]
  condition <- intersection$Condition[i]
  pbmc_integre$condition[pbmc_integre$HTO_classification == hto] <- condition
}
metadata <- pbmc_integre@meta.data
levels(metadata$condition)
metadata$condition <- factor(metadata$condition, levels = c("HV", "MS", "Ag_MS"))
pbmc_integre@meta.data <- metadata

#Création d'une colonne pour trier par run
pbmc_integre$run <- NA
pbmc_integre$run <- substr(pbmc_integre$orig.ident, 1, 4)

#Création d'une colonne pour trier par patient
pbmc_integre$patient <- paste(pbmc_integre@meta.data$orig.ident, pbmc_integre$condition, sep = "_")

#Création d'une colonne pour le sexe des patients
pbmc_integre@meta.data <-  cbind(pbmc_integre@meta.data, "Sexe"=mapvalues(x=pbmc_integre@meta.data$patient, 
                                                                    from=c("run1p1_MS", "run1p1_HV", "run1p1_Ag_MS", "run1p2_HV", "run1p2_Ag_MS", "run1p3_MS", "run1p3_HV", "run2p1_HV", "run2p1_Ag_MS", "run2p1_MS", "run2p2_HV", "run2p2_MS", "run3p1_HV", "run3p1_MS", "run3p1_Ag_MS", "run4p1_HV", "run4p1_MS", "run4p1_Ag_MS", "run3p2_MS", "run2p3_MS", "run2p3_HV", "run4p2_Ag_MS", "run4p3_Ag_MS", "run1p2_MS", "run1p3_Ag_MS", "run2p2_Ag_MS", "run3p2_Ag_MS", "run3p2_HV", "run2p3_Ag_MS", "run4p2_MS", "run4p2_HV", "run4p3_MS", "run4p3_HV"), 
                                                                    to=c("F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H")))
#Observation de la UMAP annotée
Idents(pbmc_integre) <- "populations"
DimPlot(object = pbmc_integre, reduction = "umap", label = TRUE, label.size = 7)

#Comparaison de la fréquence des cellules dans chaque sous clusters
df <- as.data.frame(table(pbmc_integre@meta.data[,c("condition", "populations")]))
ggplot(df, aes(x=condition, y=Freq, fill = populations))+geom_bar(stat="identity", position ="fill")
contingency_table <- table(pbmc_integre@meta.data$condition, pbmc_integre@meta.data$populations)
chisq.test(contingency_table)

#Subset des monocytes
mono <- subset(pbmc_integre, populations %in% "Mono")
saveRDS(mono, file = "monocyte.RDS")

#Chargement du dataset contenant uniquement les monocytes
mono <- readRDS("monocyte.RDS")
DefaultAssay(mono) <- "integrated"

#Observation des analyses primaires déjà effectué
VlnPlot(mono, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =3, pt.size = 0)
plot1 <- FeatureScatter(mono, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mono, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#PCA
mono <- ScaleData(mono, verbose = FALSE)
mono <- FindVariableFeatures(mono, selection.method = "vst", nfeatures = 4000)  
mono <- RunPCA(mono, npcs = 30, verbose = TRUE)
DimPlot(mono, reduction = "pca", group.by = "run")
DimPlot(mono, reduction = "pca", group.by = "condition")
DimHeatmap(mono, dims = 1:9, cells = 500, balanced = TRUE)

#Elbowplot
ElbowPlot(mono)

#UMAP
mono <- FindNeighbors(mono, dims = 1:20)
mono <- FindClusters(mono, resolution = 1)
mono <- RunUMAP(mono, dims = 1:20, n.neighbors = 50, min.dist = 0.2, n.components = 2)
UMAPPlot(mono, reduction = "umap")

#Détermination des marqueurs
Idents(mono) <- "seurat_clusters"
mono.markers <- FindAllMarkers(mono, only.pos = TRUE, min.pct = 0.25)

#Visualisation des Topmarkers sous forme de dotplot
genes_signif <- subset(mono.markers, mono.markers$p_val_adj < 0.05)
top_markers <- genes_signif %>%
  top_n(50, avg_log2FC)
dittoDotPlot(mono,vars=unique(top_markers$gene), group.by= "seurat_clusters") +  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

#Visualisation des marqueurs d'intérêts sous forme de dotplot
features <- c( "FCGR1A","STAT1", "IFIT1", "ISG15", "HERC5", "CD14","LYZ", "FCN1", "S100A8", "S100A9", "S100A12", "CCL5", "SDPR", "TUBB1", "CD79B", "TCF7L2","VNN2","SELL", "ITGAM", "FCGR3A","MS4A7","RHOC","CTSL", "FCER1A", "CLEC10A", "CD1C", "CD19", "MS4A1", "CD79A", "CD40", "PAX5","IL4R", "NCAM1", "CD8A", "CD8B","CD27", "KLRG1", "NKG7","GNLY", "CD247","GZMA", "GZMB", "IL7R", "CCR7")
DotPlot(mono, features = features) +  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + RotatedAxis()

#Visualisation des marqueurs d'intérêts sous forme de featureplot
FeaturePlot(mono, reduction = "umap", features=c("CD14","LYZ", "S100A8", "SELL"), min.cutoff = 'q10', label = TRUE) #Monocyte classique
FeaturePlot(mono, reduction = "umap", features=c("STAT1", "IFIT1", "ISG15", "HERC5"), min.cutoff = 'q10', label = TRUE) #MC-GSI
FeaturePlot(mono, reduction = "umap", features=c("FCGR3A", "RHOC"), min.cutoff = 'q10', label = TRUE) #Monocyte non classiques
FeaturePlot(mono, reduction = "umap", features=c("FCER1A","CD1C", "CLEC10A"), min.cutoff = 'q10', label = TRUE) #DC
FeaturePlot(mono, reduction = "umap", features=c("CD19", "MS4A1", "CD79A", "PAX5"), min.cutoff = 'q10', label = TRUE) #LB
FeaturePlot(mono, reduction = "umap", features=c("NCAM1", "KLRG1", "NKG7", "IL7R"), min.cutoff = 'q10', label = TRUE) #NKT

#Visualisation des marqueurs d'intérêts sous forme de VlnPLot
DefaultAssay(mono) <- "RNA"
VlnPlot(mono, features=c("IFIT1", "ISG15", "HERC5"), pt.size = 0) #MC-GSI
VlnPlot(mono, features=c("CD19", "MS4A1", "CD79A"), pt.size = 0) #LB
VlnPlot(mono, features=c("CD8A", "NKG7", "CD3E"), pt.size = 0) #NKT
VlnPlot(mono, features=c("FCER1A","CD1C", "CLEC10A"), pt.size = 0) #cDC
VlnPlot(mono, features=c("FCGR3A","CD14"), pt.size = 0) #Monocyte

#ScatterPlot des marqueurs CD14 et CD16
dittoScatterPlot(mono, 
                 x.var = "CD14", 
                 y.var="FCGR3A",  
                 scale_color_gradient2(low = "green", mid = "gray", high = "red", midpoint = 0),
                 assay.x = "RNA", 
                 assay.y = "RNA", 
                 color.var = "seurat_clusters") +
  ggtitle("Scatterplot des marqueurs CD14 et CD16 en fonction des clusters")

#Annotation des clusters
mono@meta.data <-  cbind(mono@meta.data, "Annotation"=mapvalues(x=mono@meta.data$seurat_clusters, 
                                                                from=c(0,1,2,3,4,5,6,7,8, 9, 10, 11, 12, 13, 14, 15, 16), to=c("MC-Inf", "MC-IFN", "MC-Inf", "MC", "MC-proInf", "MNC", "MC-proInf", "MC", "MC-Inf", "LB/Mono", "MC-Inf", "LT/NK", "MI", "LT/NK", "cDC", "LB/Mono", "LB")))

#Elimination des cellules qui ne sont pas des monocytes
monocytes <- subset(mono, Annotation %in% c("MC", "MC-IFN", "MC-Inf", "MC-proInf", "MNC", "MI"))

#Enregistrement du cluster de monocytes
saveRDS(monocytes, file = "monocytes.rds")

#Chargement de l'élément seurat contenant les monocytes
setwd("/LAB-DATA/BiRD/shares/CRTI/EQ5/SC_pbmc/M1_amandine/")
monocytes <- readRDS("monocytes.rds")
DefaultAssay(monocytes) <- "integrated"

#PCA et UMAP des monocytes triés
monocytes <- RunPCA(monocytes, npcs = 30, verbose = TRUE)
monocytes <- FindNeighbors(monocytes, dims = 1:20)
monocytes <- FindClusters(monocytes, resolution = 0.4)
monocytes <- RunUMAP(monocytes, dims = 1:20, n.neighbors = 50, min.dist = 0.2, n.components = 2)
UMAPPlot(monocytes, reduction = "umap", label = TRUE)

#Visualisation des marqueurs d'intérêts sous forme de dotplot
features <- c( "FCGR1A","STAT1", "IFIT1", "ISG15", "HERC5", "CD14","LYZ", "FCN1", "S100A8", "S100A9", "S100A12", "CCL5", "SDPR", "TUBB1", "CD79B", "TCF7L2","VNN2","SELL", "ITGAM", "FCGR3A","MS4A7","RHOC","CTSL", "FCER1A", "CLEC10A", "CD1C", "CD19", "MS4A1", "CD79A", "CD40", "PAX5","IL4R", "NCAM1", "CD8A", "CD8B","CD27", "KLRG1", "NKG7","GNLY", "CD247","GZMA", "GZMB", "IL7R", "CCR7")
DotPlot(monocytes, features = features) +  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + RotatedAxis()

#Annotation des clusters
monocytes@meta.data <-  cbind(monocytes@meta.data, "Annotation2"=mapvalues(x=monocytes@meta.data$seurat_clusters, 
                                                                from=c(0,1,2,3,4,5,6), to=c("MC-proInf", "MC-proInf", "MC-Inf", "MC-Inf", "MC-GSI", "MNC", "MI")))

#Visualisation avec annotation
Idents(monocytes) <- "Annotation2"
UMAPPlot(monocytes, reduction = "umap", group.by = "Annotation2", label = TRUE, label.size = 7)
UMAPPlot(monocytes, reduction = "umap", group.by = "Annotation2", split.by = "condition", label = TRUE, label.size = 5)
features <- c( "FCGR1A","STAT1", "IFIT1", "ISG15", "HERC5", "CD14","LYZ", "FCN1", "S100A8", "S100A9", "S100A12","VNN2","SELL", "ITGAM", "FCGR3A","MS4A7","RHOC","CTSL", "CD79B", "TCF7L2")
DotPlot(monocytes, features = features) +  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + RotatedAxis()
DefaultAssay(monocytes) <- "RNA"
VlnPlot(monocytes, features=c("IFIT1", "ISG15", "HERC5"), pt.size = 0)
VlnPlot(monocytes, features=c("FCGR3A","CD14"), pt.size = 0)
VlnPlot(monocytes, features=c("S100A12","VNN2", "SELL"), pt.size = 0)

#Comparaison du nombre de cellules dans chaque sous clusters
df <- as.data.frame(table(monocytes@meta.data[,c("condition", "Annotation2")]))
ggplot(df, aes(x=condition, y=Freq, fill = Annotation2))+geom_bar(stat="identity", position ="fill")
contingency_table <- table(monocytes@meta.data$condition, monocytes@meta.data$Annotation2)
khi_test <- chisq.test(contingency_table)

#Comparaison de la fréquence de chaquesous-clusters de monocytes en fonction des conditions et par patient
df_mono <- as.data.frame(table(pbmc_integre@meta.data[,c("patient", "populations")]))
df_mono <- subset(df_mono, populations %in% "Mono")
df <- as.data.frame(table(monocytes@meta.data[,c("patient", "Annotation2")]))
merged_data <- df %>%
  left_join(df_mono %>% select(patient, Freq), by = "patient", suffix = c("", "_total"))
merged_data <- merged_data %>%
  mutate(Percentage = (Freq / Freq_total) * 100)
merged_data$condition <- substring(merged_data$patient, 8)
levels(merged_data$condition)
merged_data$condition <- factor(merged_data$condition, levels = c("HV", "MS", "Ag_MS"))

#Visualisation sous forme de Boxplot
ggplot(merged_data, aes(x = condition, y = Percentage)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, aes(color = condition)) +
  scale_color_manual(values=c("green4", "dodgerblue1","mediumpurple2")) +
  labs(title = "Fréquences de chaque sous-cluster de monocytes en fonction des conditions",
       x = "Condition",
       y = "Fréquence") +
  theme_minimal()+
  facet_wrap(~ Annotation2)

#Test statistique de Wilcoxon pour vérifier la significativité des différences de fréquence entre les conditions
perform_wilcox_test <- function(data, group_col, value_col) {
  pairwise.wilcox.test(data[[value_col]], data[[group_col]], p.adjust.method = "bonferroni")
}
results <- list()
for (cell_type in unique(merged_data$Annotation2)) {
  subset_data <- subset(merged_data, Annotation2 == cell_type)
  test_result <- perform_wilcox_test(subset_data, "condition", "Percentage")
  results[[cell_type]] <- test_result
}
results

#Fonctions pour le volcano plot
get_df <- function(results,n_genes, fc, pval, ptype){
  df = data.frame(results)
  df <- df[order(df[[ptype]]),]
  df$diffexp <- NA
  df$diffexp[df$avg_log2FC > fc & df[[ptype]] < pval] <- "UP"
  df$diffexp[df$avg_log2FC < -fc & df[[ptype]] < pval] <- "DOWN"
  df$topgenes <- NA
  deg_id <- which(!is.na(df$diffexp))
  if (length(deg_id) < n_genes) {
    df$topgenes[deg_id] <- row.names(df)[deg_id]
  } else {
    df$topgenes[deg_id[1:n_genes]] <- rownames(df)[deg_id[1:n_genes]]
  }
  return(df)
}

generate_significant_genes_table <- function(comparison_df) {
  significant_genes <- comparison_df[!is.na(comparison_df$diffexp), ]
  significant_genes <- significant_genes[order(significant_genes$p_val_adj), ]
  significant_genes_table <- data.frame(
    Gene = rownames(significant_genes),
    Average_log2FC = significant_genes$avg_log2FC,
    Adjusted_p_value = significant_genes$p_val_adj,
    Differential_expression = significant_genes$diffexp)
  return(significant_genes_table)
}

generate_volcano_plot_and_gene_table <- function(pbmc_data, condition1, condition2, log2FC_threshold, p_val_adj_threshold, title, save, save2, x) {
  Idents(pbmc_data) <- "condition"
  markers <- FindMarkers(pbmc_data, ident.1 = condition1, ident.2 = condition2, min.pct = 0.05, logfc.threshold = 0.25)
  filtered_markers <- get_df(markers, 30, log2FC_threshold, p_val_adj_threshold, 'p_val_adj')
  filtered_markers <- filtered_markers[!is.na(filtered_markers$p_val_adj), ]
  vlcoP <- ggplot(data = filtered_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexp)) + 
    geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(p_val_adj_threshold), col = "gray", linetype = 'dashed') +
    ggtitle(title) + 
    geom_point(size = 1) + theme_minimal() + 
    theme(plot.title = element_text(size = 8, hjust = 0.5)) + 
    labs(col = "Différence d'expression des gènes", x = expression("Log"[2] * "FC"), y = expression("-log10 (pValue adjusted)")) +
    scale_color_manual(values = c("#00AFBB", "#bb0c00"), labels = c("Downregulated", "Upregulated")) +
    coord_cartesian(ylim = c(0, x), xlim = c(-2.5, 2.5)) +
    theme(plot.title = element_text(hjust = 0.5, size = 15)) +
    geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = topgenes), max.overlaps = Inf)
  show(vlcoP)
  ggsave(save, vlcoP)
  significant_genes_table <- generate_significant_genes_table(filtered_markers)
  write.csv(significant_genes_table, save2)
  return(significant_genes_table)
}

# Volcano plot pseudobulk des monocytes
DefaultAssay(monocytes) <- "RNA"
mono_Ag_MS_vs_MS <- generate_volcano_plot_and_gene_table(monocytes, "Ag_MS", "MS", 0.5, 0.05, 'Volcano Plot des monocytes des patients AgMS vs MS', 'VP_mono_AgM_MS.png', 'Gene_mono_AgM_MS.csv', 90)
mono_Ag_MS_vs_HV <- generate_volcano_plot_and_gene_table(monocytes, "Ag_MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes des patients AgMS vs HV', 'VP_mono_AgM_HV.png', 'Gene_mono_AgM_HV.csv', 125)
mono_All_MS_vs_HV <- generate_volcano_plot_and_gene_table(monocytes, c("Ag_MS", "MS"), "HV", 0.5, 0.05, 'Volcano Plot des monocytes des patients All MS vs HV', 'VP_mono_All_MS_HV.png', 'Gene_mono_All_MS_HV.csv', 70)
mono_MS_vs_HV <- generate_volcano_plot_and_gene_table(monocytes, "MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes des patients MS vs HV', 'VP_mono_MS_HV.png', 'Gene_mono_MS_HV.csv', 150)

# Volcano plot de tout les monocytes classiques
all_MC <- subset(monocytes, Annotation2 %in% c("MC-proInf","MC", "MC-IFN"))
all_MC_Ag_MS_vs_MS <- generate_volcano_plot_and_gene_table(all_MC, "Ag_MS", "MS", 0.5, 0.05, 'Volcano Plot de tout les monocytes classiques des patients AgMS vs MS', 'VP_all_MC_AgMS_MS.png', 'Gene_all_MC_AgMS_MS.csv', 100)
all_MC_Ag_MS_vs_HV <- generate_volcano_plot_and_gene_table(all_MC, "Ag_MS", "HV", 0.5, 0.05, 'Volcano Plot de tout les monocytes classiques des patients AgMS vs HV', 'VP_all_MC_AgMS_HV.png', 'Gene_all_MC_AgMS_HV.csv', 120)
all_MC_All_MS_vs_HV <- generate_volcano_plot_and_gene_table(all_MC, c("Ag_MS", "MS"), "HV", 0.6, 0.05, 'Volcano Plot de tout les monocytes classiques des patients All MS vs HV', 'VP_all_MC_allMS_HV.png', 'Gene_all_MC_allMS_HV.csv', 100)
all_MC_MS_vs_HV <- generate_volcano_plot_and_gene_table(all_MC, "MS", "HV", 0.5, 0.05, 'Volcano Plot de tout les monocytes classiques des patients MS vs HV', 'VP_all_MC_MS_HV.png', 'Gene_all_MC_MS_HV.csv', 100)

# Volcano plot des monocytes classiques non IFN
allMC <- subset(monocytes, Annotation2 %in% c("MC-proInf","MC"))
allMC_proInf_Ag_MS_vs_MS <- generate_volcano_plot_and_gene_table(allMC, "Ag_MS", "MS", 0.5, 0.05, 'Volcano Plot des monocytes classiques et pro-inflammatoires des patients AgMS vs MS', 'VP_allMC_AgMS_MS.png', 'Gene_allMC_AgMS_MS.csv', 100)
allMC_proInf_Ag_MS_vs_HV <- generate_volcano_plot_and_gene_table(allMC, "Ag_MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes classiques des patients AgMS vs HV', 'VP_allMC_AgMS_HV.png', 'Gene_allMC_AgMS_HV.csv', 120)
allMC_proInf_All_MS_vs_HV <- generate_volcano_plot_and_gene_table(allMC, c("Ag_MS", "MS"), "HV", 0.5, 0.05, 'Volcano Plot des monocytes classiques et pro-inflammatoires des patients All MS vs HV', 'VP_allMC_allMS_HV.png', 'Gene_allMC_allMS_HV.csv', 65)
allMC_proInf_MS_vs_HV <- generate_volcano_plot_and_gene_table(allMC, "MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes classiques des patients MS vs HV', 'VP_all_MC_MS_HV.png', 'Gene_allMC_MS_HV.csv', 100)

# Volcano plot des monocytes pro-inflammatoires
MC_proInf <- subset(monocytes, Annotation2 %in% "MC-proInf")
MC_proInf_Ag_MS_vs_MS <- generate_volcano_plot_and_gene_table(MC_proInf, "Ag_MS", "MS", 0.5, 0.05, 'Volcano Plot des monocytes pro-inflammatoires des patients AgMS vs MS', 'VP_mono_proInf_AgMS_MS.png', 'Gene_mono_proInf_AgMS_MS.csv', 45)
MC_proInf_Ag_MS_vs_HV <- generate_volcano_plot_and_gene_table(MC_proInf, "Ag_MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes pro-inflammatoires des patients AgMS vs HV', 'VP_mono_proInf_AgMS_HV.png', 'Gene_mono_proInf_AgMS_HV.csv', 55)
MC_proInf_All_MS_vs_HV <- generate_volcano_plot_and_gene_table(MC_proInf, c("Ag_MS", "MS"), "HV", 0.5, 0.05, 'Volcano Plot des monocytes pro-inflammatoires des patients All MS vs HV', 'VP_mono_proInf_allMS_HV.png', 'Gene_mono_proInf_allMS_HV.csv', 60)
MC_proInf_MS_vs_HV <- generate_volcano_plot_and_gene_table(MC_proInf, "MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes pro-inflammatoires des patients MS vs HV', 'VP_mono_proInf_MS_HV.png', 'Gene_mono_proInf_MS_HV.csv', 50)

# Volcano plot des monocytes classiques
MC <- subset(monocytes, Annotation2 %in%  "MC")
MC_Ag_MS_vs_MS <- generate_volcano_plot_and_gene_table(MC, "Ag_MS", "MS", 0.5, 0.05, 'Volcano Plot des monocytes classiques des patients AgMS vs MS', 'VP_mono_classiques_AgMS_MS.png', 'Gene_mono_classiques_AgMS_MS.csv', 35)
MC_Ag_MS_vs_HV <- generate_volcano_plot_and_gene_table(MC, "Ag_MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes classiques des patients AgMS vs HV', 'VP_mono_classiques_AgMS_HV.png', 'Gene_mono_classiques_AgMS_HV.csv', 33)
MC_All_MS_vs_HV <- generate_volcano_plot_and_gene_table(MC, c("Ag_MS", "MS"), "HV", 0.5, 0.05, 'Volcano Plot des monocytes classiques des patients All MS vs HV', 'VP_mono_classiques_AllMS_HV.png', 'Gene_mono_classiques_AllMS_HV.csv', 35)
MC_MS_vs_HV <- generate_volcano_plot_and_gene_table(MC, "MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes classiques des patients MS vs HV', 'VP_mono_classique_MS_HV.png', 'Gene_mono_classiques_MS_HV.csv', 35)

# Volcano plot des monocytes intermédiaires
MI <- subset(monocytes, Annotation %in% "MI")
MI_Ag_MS_vs_MS <- generate_volcano_plot_and_gene_table(MI, "Ag_MS", "MS", 0.5, 0.05, 'Volcano Plot des monocytes intermédiaires des patients AgMS vs MS', 'VP_mono_intermédiaires_AgMS_MS.png', 'Gene_mono_intermédiaires_AgMS_MS.csv', 5)
MI_Ag_MS_vs_HV <- generate_volcano_plot_and_gene_table(MI, "Ag_MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes intermédiaires des patients AgMS vs HV', 'VP_mono_intermédiaires_AgMS_HV.png', 'Gene_mono_intermédiaires_AgMS_HV.csv', 10)
MI_All_MS_vs_HV <- generate_volcano_plot_and_gene_table(MI, c("Ag_MS", "MS"), "HV", 0.5, 0.05, 'Volcano Plot des monocytes intermédiaires des patients All MS vs HV', 'VP_mono_intermédiaires_AllMS_HV.png', 'Gene_mono_intermédiaires_AllMS_HV.csv', 10)
MI_MS_vs_HV <- generate_volcano_plot_and_gene_table(MI, "MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes intermédiaires des patients MS vs HV', 'VP_mono_intermédiaires_MS_HV.png', 'Gene_mono_intermédiaires_MS_HV.csv', 7)

# Volcano plot des monocytes non classiques
MNC <- subset(monocytes, Annotation == "MNC")
MNC_Ag_MS_vs_MS <- generate_volcano_plot_and_gene_table(MNC, "Ag_MS", "MS", 0.5, 0.05, 'Volcano Plot des monocytes non classiques des patients AgMS vs MS', 'VP_mono_NC_AgMS_MS.png', 'Gene_mono_NC_AgMS_MS.csv', 13)
MNC_Ag_MS_vs_HV <- generate_volcano_plot_and_gene_table(MNC, "Ag_MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes non classiques des patients AgMS vs HV', 'VP_mono_NC_AgMS_HV.png', 'Gene_mono_NC_AgMS_HV.csv', 17)
MNC_All_MS_vs_HV <- generate_volcano_plot_and_gene_table(MNC, c("Ag_MS", "MS"), "HV", 0.5, 0.05, 'Volcano Plot des monocytes non classiques des patients All MS vs HV', 'VP_mono_NC_AllMS_HV.png', 'Gene_mono_NC_AllMS_HV.csv', 10)
MNC_MS_vs_HV <- generate_volcano_plot_and_gene_table(MNC, "MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes non classiques des patients MS vs HV', 'VP_mono_NC_MS_HV.png', 'Gene_mono_NC_MS_HV.csv', 15)

# Volcano plot des monocytes voie IFN
MC_IFN <- subset(monocytes, Annotation == "MC-IFN")
MC_IFN_Ag_MS_vs_MS <- generate_volcano_plot_and_gene_table(MC_IFN, "Ag_MS", "MS", 0.5, 0.05, 'Volcano Plot des monocytes classiques IFN des patients AgMS vs MS', 'VP_MC_IFN_AgMS_MS.png', 'Gene_MC_IFN_AgMS_MS.csv', 10)
MC_IFN_Ag_MS_vs_HV <- generate_volcano_plot_and_gene_table(MC_IFN, "Ag_MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes classiques IFN des patients AgMS vs HV', 'VP_MC_IFN_AgMS_HV.png', 'Gene_MC_IFN_AgMS_HV.csv', 50)
MC_IFN_All_MS_vs_HV <- generate_volcano_plot_and_gene_table(MC_IFN, c("Ag_MS", "MS"), "HV", 0.5, 0.05, 'Volcano Plot des monocytes classiques IFN des patients All MS vs HV', 'VP_MC_IFN_AllMS_HV.png', 'Gene_MC_IFN_AllMS_HV.csv', 45)
MC_IFN_MS_vs_HV <- generate_volcano_plot_and_gene_table(MC_IFN, "MS", "HV", 0.5, 0.05, 'Volcano Plot des monocytes classiques IFN des patients MS vs HV', 'VP_MC_IFN_MS_HV.png', 'Gene_MC_IFN_MS_HV.csv', 45)

#Enrichissement via EnrichR
dbs_of_interest <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022")
enrichment_results <- enrichr(gene_allMC_AgMS_MS$Gene, dbs_of_interest)
enrichment_results <- enrichr(gene_allMC_AgMS_HV$Gene, dbs_of_interest)
enrichment_results <- enrichr(gene_allMC_AllMS_HV$Gene, dbs_of_interest)
enrichment_results <- enrichr(gene_allMC_MS_HV$Gene, dbs_of_interest)
plotEnrich(enrichment_results[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enrichment_results[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

#VlnPlot de l'expression de marqueurs d'intérêts
VlnPlot(monocytes, features = "TSC22D3", group.by="Annotation2", split.by = "condition", pt.size=0, log = TRUE, cols=c("green4", "dodgerblue1","mediumpurple2"))
VlnPlot(monocytes, features = "HDAC9", group.by="Annotation2", split.by = "condition", pt.size=0, log = TRUE, cols=c("green4", "dodgerblue1","mediumpurple2"))
VlnPlot(monocytes, features = "FKBP5", group.by="Annotation2", split.by = "condition", pt.size=0, log = TRUE, cols=c("green4", "dodgerblue1","mediumpurple2"))
VlnPlot(monocytes, features = "KLF9", group.by="Annotation2", split.by = "condition", pt.size=0, log = TRUE, cols=c("green4", "dodgerblue1","mediumpurple2"))
VlnPlot(monocytes, features = "SAP30", group.by="Annotation2", split.by = "condition", pt.size=0, log = TRUE, cols=c("green4", "dodgerblue1","mediumpurple2"))
VlnPlot(monocytes, features = "PDK4", group.by="Annotation2", split.by = "condition", pt.size=0, log = TRUE, cols=c("green4", "dodgerblue1","mediumpurple2"))
VlnPlot(monocytes, features = "CXCR4", group.by="Annotation2", split.by = "condition", pt.size=0, log = TRUE, cols=c("green4", "dodgerblue1","mediumpurple2"))

#WGCNA
enableWGCNAThreads(nThreads = 8)
DefaultAssay(monocytes) <- "RNA"
monocytes <- RunHarmony(monocytes, group.by.vars = "Annotation")
mono_wgcna <- SetupForWGCNA(monocytes, gene_select="fraction", fraction = 0.05, wgcna_name = "monocytes")
mono_wgcna <- MetacellsByGroups(seurat_obj = mono_wgcna, group.by = c("Annotation", "condition"), reduction = 'harmony', k = 25, max_shared = 10, ident.group = 'Annotation')
mono_wgcna <- NormalizeMetacells(mono_wgcna)
MC_proInf_wgcna <- SetDatExpr(mono_wgcna, group_name = c("MC-proInf","MC","MC-Inf"), group.by="Annotation", assay='RNA', slot='data')
MC_proInf_wgcna <- TestSoftPowers(MC_proInf_wgcna, networkType = 'signed')
plot_list <- PlotSoftPowers(MC_proInf_wgcna)
wrap_plots(plot_list, ncol=2)
MC_proInf_wgcna <- ConstructNetwork(MC_proInf_wgcna, tom_name="all-MC", overwrite_tom = TRUE)
PlotDendrogram(MC_proInf_wgcna, main='INH hdWGCNA Dendrogram')
MC_proInf_wgcna <- ModuleEigengenes(MC_proInf_wgcna, group.by.vars="condition")
MEs <- MC_proInf_wgcna$eigengenes 
DefaultAssay(MC_proInf_wgcna) <- "RNA"
hMEs <- GetMEs(MC_proInf_wgcna)
MC_proInf_wgcna <- ModuleConnectivity(MC_proInf_wgcna, group.by = 'Annotation', group_names = c("MC-proInf","MC","MC-Inf"))
MC_proInf_wgcna <- ResetModuleNames(MC_proInf_wgcna, new_name = "MC-Mod")
PlotKMEs(MC_proInf_wgcna, ncol=5)

#Etude différencié homme:femme
monocytes@meta.data <-  cbind(monocytes@meta.data, "Sexe"=mapvalues(x=monocytes@meta.data$patient, 
                                                                           from=c("run1p1_MS", "run1p1_HV", "run1p1_Ag_MS", "run1p2_HV", "run1p2_Ag_MS", "run1p3_MS", "run1p3_HV", "run2p1_HV", "run2p1_Ag_MS", "run2p1_MS", "run2p2_HV", "run2p2_MS", "run3p1_HV", "run3p1_MS", "run3p1_Ag_MS", "run4p1_HV", "run4p1_MS", "run4p1_Ag_MS", "run3p2_MS", "run2p3_MS", "run2p3_HV", "run4p2_Ag_MS", "run4p3_Ag_MS", "run1p2_MS", "run1p3_Ag_MS", "run2p2_Ag_MS", "run3p2_Ag_MS", "run3p2_HV", "run2p3_Ag_MS", "run4p2_MS", "run4p2_HV", "run4p3_MS", "run4p3_HV"), 
to=c("F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H")))

df <- as.data.frame(table(monocytes@meta.data[,c("Sexe", "condition", "Annotation2")]))
ggplot(df, aes(x=Sexe, y=Freq, fill = Annotation2)) + 
  geom_bar(stat="identity", position ="fill") +
  facet_wrap(~ condition) +
  labs(y = "Proportion", x = "Sexe") + 
  theme_minimal()
contingency_table <- table(mono@meta.data$condition, mono@meta.data$Annotation2)
khi_test <- chisq.test(contingency_table)

#Etude via DESeq2
Idents(monocytes) <- "condition"
pseudo_monocytes <- AggregateExpression(monocytes, assays = 'RNA', slot = 'count', return.seurat = T, group.by = c("condition", "patient", "Annotation2"))
pseudo_monocytes$Annaconda <- paste(pseudo_monocytes$Annotation2, pseudo_monocytes$condition, sep = "_")
Idents(pseudo_monocytes) <- "condition"
markers <- FindMarkers(pseudo_monocytes, ident.1 = "Ag-MS", ident.2 = "HV", test.use = "DESeq2")
