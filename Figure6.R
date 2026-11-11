#Fig6
library(ggrepel)
library(dplyr)
library(tidyr)
library(Seurat)
library(DESeq2)
library(pheatmap)
library(ggplot2)

#Fig6A
FAS<-subset(ABC,subset=genotype %in% "fas")
WT<-subset(ABC,subset=genotype %in% "wt")
DKO<-subset(ABC,subset=genotype %in% "p38")

counts <- GetAssayData(FAS, slot = "counts")
meta <- FAS@meta.data
meta$cell_id <- colnames(FAS)
mouse <- unique(meta$mouse)
pb_list <- lapply(mouse, function(s) {
  cells <- meta %>% filter(mouse == s) %>% pull(cell_id)
  if (length(cells) == 1) {
    counts[, cells]
  } else {
    Matrix::rowSums(counts[, cells])
  }
})
names(pb_list) <- mouse
pb_counts <- do.call(cbind, pb_list)
colnames(pb_counts) <- mouse
FAS_pb_list<-pb_list
FAS_pb_counts<-pb_counts

counts <- GetAssayData(WT, slot = "counts")
meta <- WT@meta.data
meta$cell_id <- colnames(WT)
mouse <- unique(meta$mouse)
pb_list <- lapply(mouse, function(s) {
  cells <- meta %>% filter(mouse == s) %>% pull(cell_id)
  if (length(cells) == 1) {
    counts[, cells]
  } else {
    Matrix::rowSums(counts[, cells])
  }
})
names(pb_list) <- mouse
pb_counts <- do.call(cbind, pb_list)
colnames(pb_counts) <- mouse
WT_pb_list<-pb_list
WT_pb_counts<-pb_counts

counts <- GetAssayData(DKO, slot = "counts")
meta <- DKO@meta.data
meta$cell_id <- colnames(DKO)
mouse <- unique(meta$mouse)
pb_list <- lapply(mouse, function(s) {
  cells <- meta %>% filter(mouse == s) %>% pull(cell_id)
  if (length(cells) == 1) {
    counts[, cells]
  } else {
    Matrix::rowSums(counts[, cells])
  }
})
names(pb_list) <- mouse
pb_counts <- do.call(cbind, pb_list)
colnames(pb_counts) <- mouse
DKO_pb_list<-pb_list
DKO_pb_counts<-pb_counts
mat<-cbind(WT_pb_counts,FAS_pb_counts,DKO_pb_counts)
genes <- c(
  "Ramp3", "Gpc1", "Ndnf", "Cd70", "P2ry13",
  "Timd2", "Hspg2", "Procr", "Nbl1", "Cxcl10",
  "Spon1", "Tnfrsf8", "Dkk3", "Tlr12", "Itgax",
  "Fgl2", "Prg2", "Cd276"
)
genes_in_mat <- intersect(genes, rownames(mat))
if (length(genes_in_mat) == 0) stop("no gene！")
mat_selected <- mat[genes_in_mat, , drop = FALSE] 
mat_scaled <- t(scale(t(mat_selected))) 
col_order <- c("WT1", "WT2", "WT3", "FAS1", "FAS2", "FAS3", "DKO1", "DKO2", "DKO3")
col_order <- col_order[col_order %in% colnames(mat_scaled)]
mat_ordered <- mat_scaled[, col_order]
mat_clipped <- pmax(pmin(mat_ordered, 1.2), -1.2)
heatmap <- pheatmap(
  mat_clipped,  
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("navyblue", "white", "firebrick3"))(200),
  border_color = NA,
  fontsize = 10,
  main = "Selected Genes Expression Heatmap",
  angle_col = 0,
  fontsize_col = 10,
  color_breaks = seq(-1.2, 1.2, length.out = 201) 
)

#Fig6B
sample_names <- colnames(mat)
group <- ifelse(grepl("^FAS", sample_names), "FAS", "WT")
meta <- data.frame(
  sample = sample_names,
  group = factor(group, levels = c("FAS", "WT"))  
)
rownames(meta) <- sample_names
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = meta,
                              design = ~ group)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "FAS", "WT"))
res_df <- as.data.frame(res[order(res$padj), ])
DEGs<-res_df
DEGs$label <- ifelse(DEGs$gene %in% genes, "Highlighted", "Not Highlighted")
DEGs$ "-log10pvalue" <- -log10( DEGs$pvalue)
cut_off_pvalue =0.05 
cut_off_log2FC =1 
DEGs$Sig = ifelse(DEGs$pvalue< cut_off_pvalue &    
                    abs(DEGs$log2FoldChange) >= cut_off_log2FC, 
                  ifelse(DEGs$log2FoldChange > cut_off_log2FC ,'Up','Down'),'no')
DEGs <- DEGs %>%
  mutate(color_group = case_when(
    label == "Highlighted" & log2FoldChange > 0 ~ "Up",
    label == "Highlighted" & log2FoldChange < 0 ~ "Down",
    TRUE ~ "Not Highlighted"
  ))
p1 <- ggplot(DEGs, aes(x = log2FoldChange, y = `-log10pvalue`)) +
  geom_point(alpha = 1, size = 1.5, aes(color = color_group),
             data = subset(DEGs,color_group=="Not Highlighted")) +
  geom_point(alpha = 1, size = 3, aes(color = color_group),
             data = subset(DEGs,color_group=="Down")) +
  geom_point(alpha = 1, size = 3, aes(color = color_group),
             data = subset(DEGs,color_group=="Up"))+
  scale_color_manual(values = c("Up" = "#ff4757",     
                                "Down" = "#546de5",   
                                "Not Highlighted" = "#d2dae2"
  ),
  breaks = c("Up", "Down"),          
  guide = "legend") +               
  xlim(c(-5,5)) +
  ylim(c(0,8))+
  labs(x = "log2FC", y = "-log10Pvalue") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank())
p2 <- p1 + 
  geom_text_repel(
    data = subset(DEGs, gene %in% genes),
    aes(label = gene),
    size = 5,
    force = 1,                   
    max.iter = 1e3,               
    box.padding = 0.3,
    point.padding = 0.3,
    segment.size = 0.5,
    segment.color = "grey50",
    direction = "both",
    remove_overplot = FALSE,      
    max.overlaps = Inf,
    nudge_x = 0.1,              
    na.rm = TRUE,
    seed = 123
  )
#Fig6C
sample_names <- colnames(mat)
group <- ifelse(grepl("^FAS", sample_names), "FAS", "DKO")
meta <- data.frame(
  sample = sample_names,
  group = factor(group, levels = c("FAS", "DKO"))  
)
rownames(meta) <- sample_names
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = meta,
                              design = ~ group)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "FAS", "DKO"))
res_df <- as.data.frame(res[order(res$padj), ])
DEGs<-res_df
DEGs$label <- ifelse(DEGs$gene %in% genes, "Highlighted", "Not Highlighted")
DEGs$ "-log10pvalue" <- -log10( DEGs$pvalue)
cut_off_pvalue =0.05 
cut_off_log2FC =1 
DEGs$Sig = ifelse(DEGs$pvalue< cut_off_pvalue &    
                    abs(DEGs$log2FoldChange) >= cut_off_log2FC, 
                  ifelse(DEGs$log2FoldChange > cut_off_log2FC ,'Up','Down'),'no')
DEGs <- DEGs %>%
  mutate(color_group = case_when(
    label == "Highlighted" & log2FoldChange > 0 ~ "Up",
    label == "Highlighted" & log2FoldChange < 0 ~ "Down",
    TRUE ~ "Not Highlighted"
  ))
p1 <- ggplot(DEGs, aes(x = log2FoldChange, y = `-log10pvalue`)) +
  geom_point(alpha = 1, size = 1.5, aes(color = color_group),
             data = subset(DEGs,color_group=="Not Highlighted")) +
  geom_point(alpha = 1, size = 3, aes(color = color_group),
             data = subset(DEGs,color_group=="Down")) +
  geom_point(alpha = 1, size = 3, aes(color = color_group),
             data = subset(DEGs,color_group=="Up"))+
  scale_color_manual(values = c("Up" = "#ff4757",     
                                "Down" = "#546de5",   
                                "Not Highlighted" = "#d2dae2"
  ),
  breaks = c("Up", "Down"),          
  guide = "legend") +               
  xlim(c(-5,5)) +
  ylim(c(0,8))+
  labs(x = "log2FC", y = "-log10Pvalue") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank())
p2 <- p1 + 
  geom_text_repel(
    data = subset(DEGs, gene %in% genes),
    aes(label = gene),
    size = 5,
    force = 1,                   
    max.iter = 1e3,               
    box.padding = 0.3,
    point.padding = 0.3,
    segment.size = 0.5,
    segment.color = "grey50",
    direction = "both",
    remove_overplot = FALSE,      
    max.overlaps = Inf,
    nudge_x = 0.1,              
    na.rm = TRUE,
    seed = 123
  )
