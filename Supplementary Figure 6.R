#Supplementary Figures 6
library(Seurat)
library(Matrix)
library(DESeq2)
library(ggplot2)
library(dplyr)
#FigS6A
all_mat_seurat <- subset(all_mat_seurat, nCount_RNA > 400 & nCount_RNA < 10000 & percentage.mt< 50) # filter
dim(all_mat_seurat)
setwd("/home/user1/total/file_10x/")
write.table(data.frame(rownames(all_mat_seurat@assays$RNA@counts),rownames(all_mat_seurat@assays$RNA@counts)),file = 'genes.tsv',
            quote = F,sep = '\t',
            col.names = F,row.names = F)
write.table(colnames(all_mat_seurat@assays$RNA@counts),file = 'barcodes.tsv',quote = F,
            col.names = F,row.names = F)
counts <- as(all_mat_seurat@assays$RNA@counts, "sparseMatrix")
writeMM(counts, file = "matrix.mtx") 
total.data<- Read10X(data.dir = "/home/user1/total/file_10x/")
total <- CreateSeuratObject(counts = total.data, project = "total")
total <- NormalizeData(total, normalization.method = "LogNormalize", scale.factor = 10000)
total<- FindVariableFeatures(total, selection.method = "vst", nfeatures = 2000)
top100 <- head(VariableFeatures(total), 100)
all.genes <- rownames(total)
total <- ScaleData(total, features = all.genes)
total <- RunPCA(total, features = VariableFeatures(object = total))
DimPlot(total, reduction = "pca") + NoLegend()
print(total[["pca"]], dims = 1:10, nfeatures = 5)
DimHeatmap(total, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(total,ndims = 100)
total<- FindNeighbors(total, dims = 1:30)
total<- FindClusters(total, resolution = 0.5)
total<- RunUMAP(total, dims = 1:30)
DimPlot(total, reduction = "umap")

#FigS6B
feature_genes<-c("Cxcl10","Cxcr3")
output_dir <- "D:/FAS/FigS6/"
dir.create(output_dir, showWarnings = FALSE)
for (gene in feature_genes) {
  if (!gene %in% rownames(Fas)) {
    warning(paste("Gene", gene, "not found. Skipping..."))
    next
  }
  
  p <- FeaturePlot(Fas, features = gene, order = TRUE, pt.size = 1.2) +
    theme(
      axis.title.x = element_blank(),       
      axis.title.y = element_blank(),        
      axis.line = element_blank(),           
      axis.ticks = element_blank(),          
      axis.text = element_blank(),           
      legend.position = "none",             
      plot.title = element_blank(),          
      panel.background = element_blank(),   
      panel.grid = element_blank()#,         
      
      #text = element_text(family = "Times New Roman")  
    )
  
  print(p)
  
  output_file <- file.path(output_dir, paste0(gene, "表达水平.svg"))
  
  
  ggsave(
    filename = output_file,
    ,plot = p               
    ,device = "svg"        
    ,width = 8            
    ,height = 8             
    ,units = "in"           
    ,bg = "transparent"      
  )
  cat("Saved:", output_file, "\n")
}
#FigS6C
ABC<-subset(total,subset=annotation %in% "ABC")
Mon<-subset(total,subset=annotation %in% "Mon-Mac")
ABC$mouse <- ifelse(
  grepl("^FAS", ABC$mouse),
  paste0("fas", sub("^FAS", "", ABC$mouse)),
  paste0("p38", sub("^DKO", "", ABC$mouse))
)
Mon$mouse <- ifelse(
  grepl("^FAS", Mon$mouse),
  paste0("fas", sub("^FAS", "", Mon$mouse)),
  paste0("p38", sub("^DKO", "", Mon$mouse))
)
pb_ABC <- AggregateExpression(ABC, assays = "RNA", slot = "counts", 
                              group.by = "mouse", return.seurat = FALSE)
pb_Mon <- AggregateExpression(Mon, assays = "RNA", slot = "counts", 
                              group.by = "mouse", return.seurat = FALSE)
counts_ABC <- pb_ABC$RNA
counts_Mon <- pb_Mon$RNA
colnames(counts_ABC) <- paste0(colnames(counts_ABC), "_ABC")
colnames(counts_Mon) <- paste0(colnames(counts_Mon), "_Mon")
common_genes <- intersect(rownames(counts_ABC), rownames(counts_Mon))
counts_all <- cbind(counts_ABC[common_genes, ], counts_Mon[common_genes, ])
coldata <- data.frame(colname = colnames(counts_all), stringsAsFactors = FALSE)
coldata$mouse_id <- sub("_(ABC|Mon)$", "", coldata$colname)
coldata$celltype <- sub("^.*_(ABC|Mon)$", "\\1", coldata$colname)
coldata$genotype <- ifelse(grepl("^fas", coldata$mouse_id), "fas", "p38")
coldata$genotype <- factor(coldata$genotype, levels = c("fas", "p38"))
coldata$celltype <- factor(coldata$celltype, levels = c("Mon", "ABC"))
rownames(coldata) <- coldata$colname
coldata <- coldata[colnames(counts_all), ]
keep <- rowSums(counts_all >= 10) >= 2
counts_all <- counts_all[keep, ]
dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData   = coldata,
  design    = ~ genotype + celltype + genotype:celltype
)
dds <- DESeq(dds)
res_int <- results(dds, name = "genotypep38.celltypeABC")
res_Mon <- results(dds, contrast = c("genotype", "p38", "fas"))
res_ABC <- results(dds, contrast = list(c("genotype_p38_vs_fas", "genotypep38.celltypeABC")))
cxcl10_int <- res_int["Cxcl10", ]
cxcl10_Mon <- res_Mon["Cxcl10", ]
cxcl10_ABC <- res_ABC["Cxcl10", ]
cxcl10_table <- data.frame(
  comparison = c("Mon: p38 vs fas", "ABC: p38 vs fas", "Interaction (ABC - Mon)"),
  log2FC = c(cxcl10_Mon$log2FoldChange, cxcl10_ABC$log2FoldChange, cxcl10_int$log2FoldChange),
  lfcSE = c(cxcl10_Mon$lfcSE, cxcl10_ABC$lfcSE, cxcl10_int$lfcSE),
  pvalue = c(cxcl10_Mon$pvalue, cxcl10_ABC$pvalue, cxcl10_int$pvalue)
)
cxcl10_table$CI_low  <- cxcl10_table$log2FC - 1.96 * cxcl10_table$lfcSE
cxcl10_table$CI_high <- cxcl10_table$log2FC + 1.96 * cxcl10_table$lfcSE
cxcl10_table$sig <- ifelse(cxcl10_table$pvalue < 0.001, "***",
                           ifelse(cxcl10_table$pvalue < 0.01,  "**",
                                  ifelse(cxcl10_table$pvalue < 0.05,  "*", "ns")))
cxcl10_forest_plot <- cxcl10_table %>%
  filter(comparison != "Interaction (ABC - Mon)") %>%
  droplevels()
cxcl10_forest_plot$comparison <- factor(
  cxcl10_forest_plot$comparison,
  levels = c("ABC: p38 vs fas", "Mon: p38 vs fas")
)
p_forest <- ggplot(cxcl10_forest_plot, aes(x = log2FC, y = comparison)) +
  geom_segment(aes(x = CI_low, xend = CI_high, y = comparison, yend = comparison),
               size = 6, color = "grey85", alpha = 0.6) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), 
                 height = 0.15, color = "black") +
  geom_point(aes(color = comparison), size = 6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text(aes(x = CI_high, label = sig), 
            hjust = -0.3, vjust = 0.5, size = 6) +
  scale_color_manual(values = c(
    "ABC: p38 vs fas" = "#E64B35", 
    "Mon: p38 vs fas" = "#3C5488" 
  )) +
  scale_x_continuous(
    limits = c(-3.0, 0.8), 
    breaks = seq(-3.0, 0.5, 0.5),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  scale_y_discrete(
    labels = c("ABC: p38 vs fas" = "ABC", "Mon: p38 vs fas" = "Mon-Mac")
  ) +
  labs(
    x = "Log2FoldChange (DKO vs KO)", 
    y = NULL,
    title = "Cxcl10: DKO vs KO",
    subtitle = "Error bars: 95% CI"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank()
  )
print(p_forest)