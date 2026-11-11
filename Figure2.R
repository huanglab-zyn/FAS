#Fig2
library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)
library(svglite)
library(metacell)
library(stringr)
library(foreach)
library(pheatmap)
library(tgconfig)


custom_colors <- c(
  "B" = "#E6194B", 
  "cDC1" = "#3CB44B",
  "cDC2" = "#91A0C2",
  "Mon-Mac" = "#4363D8",
  "Mac" = "#911EB4",
  "mregDC" = "#42D4F4",
  "NK" = "#F032E6",
  "pDC" = "#BFBD00",
  "plasma" = "#fabed4",
  "LTi-like" = "#469990",
  "ABC"= "#FFA500",
  "T"="black",
  "unknown"="black",
  "#N/A" = "#FF7F7F", 
  "CD4_T" = "#FFD700", 
  "CD8_T" = "#693A60", 
  "naiveB" = "#FF69B4", 
  "DC" = "#00FA9A", 
  "proliferative T" = "#00CBA9", 
  "DN_T" = "#ADFF2F", 
  "GC B" = "#3A6940", 
  "siglecgB" = "#EE82EE",
  "Treg" = "grey50",
  "plasmaB"="red"
)
#CD11c/total
custom_colors <- c(
  "B" = "#D72348", 
  "cDC1" = "#40AB4B",
  "cDC2" = "#8E9DBD",
  "Mon-Mac" = "#405FAB",
  "Mac" = "#793690",
  "mregDC" = "#54C1DF",
  "NK" = "#AE529B",
  "pDC" = "#B8B619",
  "plasma" = "#F4BBCF",
  "LTi-like" = "#45948B",
  "ABC"= "#F4A016",
  "T"="black",
  "unknown"="black",
  "#N/A" = "#FF7F7F", 
  "CD4_T" = "#F5D01F", 
  "CD8_T" = "#66385D", 
  "naiveB" = "#E868A2", 
  "DC" = "#6ABE86",
  "proliferative T" = "#00CBA9", 
  "DN_T" = "#ACCF3E", 
  "GC B" = "#39673F", 
  "siglecgB" = "#C486B8", #
  "Treg" = "#7F8080",
  "plasmaB"="#E71F19",
  "wt"="#9DCFC8",
  "fas"="#DE9BAF",
  "p38"="#BDC7DE"
)

#Fig2B
all_mat_seurat <- subset(all_mat_seurat, nCount_RNA > 400 & nCount_RNA < 10000 & percentage.mt< 50) # filter
dim(all_mat_seurat)
setwd("/home/user1/file_10x/")
library(Matrix)
write.table(data.frame(rownames(all_mat_seurat@assays$RNA@counts),rownames(all_mat_seurat@assays$RNA@counts)),file = 'genes.tsv',
            quote = F,sep = '\t',
            col.names = F,row.names = F)
write.table(colnames(all_mat_seurat@assays$RNA@counts),file = 'barcodes.tsv',quote = F,
            col.names = F,row.names = F)

counts <- as(all_mat_seurat@assays$RNA@counts, "sparseMatrix")
writeMM(counts, file = "matrix.mtx") 


###########---------------------------------
### after python



obs<-read.csv("/home/user1/obs_CD11c-mesen.csv")
obs<-obs[obs$metacell > -1, ]

obs$metacell=obs$metacell+1


write.csv(obs,"/home/user1/obs_CD11c-mesen_filtered.csv")

keep_cells <- obs$X

igcells<- setdiff(colnames(mat@mat),keep_cells)
mcell_mat_ignore_cells(mat_id,mat_id,igcells) # remove igcells in second mat_id, and output to the first mat_id

nms <- unique(c(rownames(mat@mat), rownames(mat@ignore_gmat)))
######for mouse############
pre_nr_term<-c("^Rps", "^Rpl", "^Mrps", "^Mrpl")
pre_nr_genes <- foreach(i=pre_nr_term, .combine = c) %do% grep(i, nms, v=T)
pre_ex_genes <- c("Malat1", "Xist")
pre_bad_genes <- unique(c(pre_nr_genes, pre_ex_genes))
pre_bad_genes
print(length(pre_bad_genes))
#######################
mcell_mat_ignore_genes(new_mat_id=mat_id, mat_id=mat_id, pre_bad_genes, reverse=F)

mcmc<-as.integer(obs$metacell)
names(mcmc)<-obs$X

mat<-scdb_mat(mat_id) 

outliers<-c()
mcell_new_mc("new_mc",mc=mcmc,outliers = as.character(outliers), mat) # create a new metacell with the matrix "mat", which does not have outliers. the variable "mcmc" used to correspond "mc"
new_mc<-scdb_mc("new_mc")

new_mc

mcell_mc_add_color("new_mc",rep(NA,length(unique(mcmc))))

lfp<-round(log(new_mc@mc_fp),2)
write.csv(lfp,"/home/user1/lfp_CD11c.csv")


mcell_add_gene_stat(gstat_id="test", mat_id=mat_id, force=T)
mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="test", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="test", T_tot=100, T_top3=2) # line93-95, to find high-variable genes 

fb<-read.csv("/home/user1/mouse_fb_gene/fb_gene.csv")


gset <- scdb_gset("test_feats")
final_genes = names(gset@gene_set)
final_genes = setdiff(final_genes,fb$fb_gene) # confirming the genes used to perform clustering


gset@gene_set <- rep(1,length(final_genes))
names(gset@gene_set) <- final_genes
scdb_add_gset("test_feats",gset)
gset <- scdb_gset("test_feats")  # make the final genes list become the format that can be loaded by metacell

mcell_add_cgraph_from_mat_bknn(
  mat_id=mat_id,
  gset_id = "test_feats",
  graph_id="test_graph",
  K=40, # square root of actual cell numbers,default=40
  dsamp=T)

mcell_mgraph_logistic("new_mgraph","new_mc","test_feats")
mcell_mc2d_force_knn(mc2d_id = "new_2d",mc_id = "new_mc",graph_id = "test_graph",mgraph_id = "new_mgraph")

mc2d = scdb_mc2d("new_2d")
sc_x = as.data.frame(as.matrix(mc2d@sc_x))
nrow(sc_x)
sc_x = na.omit(sc_x)
colnames(sc_x) = c("x")
sc_y = as.data.frame(as.matrix(mc2d@sc_y))
nrow(sc_y)
sc_y = na.omit(sc_y)
colnames(sc_y) = c("y")
sc = cbind(sc_x,sc_y)
sc$Well_ID = rownames(sc)
write.csv(sc,"//home/user1/getb_sc_2d_CD11c_1020.csv")
mc_x = as.data.frame(as.matrix(mc2d@mc_x))
mc_y = as.data.frame(as.matrix(mc2d@mc_y))
mc = cbind(mc_x,mc_y)
colnames(mc) = c("mc_x","mc_y")
write.csv(mc,"//home/user1/getb_mc_2d_CD11c_1020.csv")

CD11c.data<- Read10X(data.dir = "/home/user1/file_10x/")
CD11c <- CreateSeuratObject(counts = CD11c.data, project = "CD11c")
CD11c <- NormalizeData(CD11c, normalization.method = "LogNormalize", scale.factor = 10000)
CD11c<- FindVariableFeatures(CD11c, selection.method = "vst", nfeatures = 2000)
top100 <- head(VariableFeatures(CD11c), 100)
all.genes <- rownames(CD11c)
CD11c <- ScaleData(CD11c, features = all.genes)
CD11c <- RunPCA(CD11c, features = VariableFeatures(object = CD11c))
DimPlot(CD11c, reduction = "pca") + NoLegend()
print(CD11c[["pca"]], dims = 1:10, nfeatures = 5)
DimHeatmap(CD11c, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(CD11c,ndims = 100)
CD11c<- FindNeighbors(CD11c, dims = 1:30)
CD11c<- FindClusters(CD11c, resolution = 0.5)
CD11c<- RunUMAP(CD11c, dims = 1:30)
DimPlot(CD11c, reduction = "umap")
CD11c<-subset(CD11c,subset=genotype %in% c("fas","wt"))
metadata<-CD11c@meta.data

p <- ggplot(metadata, aes(x = x, y = y, color = annotation)) +
  geom_point(size = 1.2, alpha = 1) +  
  scale_color_manual(values = custom_colors) +  
  theme_minimal() +  
  theme(
    legend.position = "none",  
    panel.grid = element_blank(),  
    axis.title = element_blank(),  
    axis.text = element_blank(),  
    axis.ticks = element_blank(),  
    plot.title = element_blank(),  
    panel.border = element_blank(),  
    plot.background = element_blank(),  
    panel.background = element_blank()  
  )
#Fig2C
feature_genes<-c("Itgax","Tbx21","Zeb2","Cd19","C1qa","Nkg7","Jchain","Xcr1","Siglech")
output_dir <- "D:/FAS/Fig2/"
dir.create(output_dir, showWarnings = FALSE)
for (gene in feature_genes) {
  if (!gene %in% rownames(CD11c)) {
    warning(paste("Gene", gene, "not found. Skipping..."))
    next
  }
  
  p <- FeaturePlot(CD11c, features = gene, order = TRUE, pt.size = 1.2) +
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

#Fig1C
p <- ggplot(metadata, aes(x = x, y = y, color = genotype)) +
  geom_point(size = 1.2, alpha = 1) +  
  scale_color_manual(values = custom_colors) +  
  theme_minimal() + 
  theme(
    legend.position = "none",  
    panel.grid = element_blank(),  
    axis.title = element_blank(),  
    axis.text = element_blank(),  
    axis.ticks = element_blank(),  
    plot.title = element_blank(),  
    panel.border = element_blank(),  
    plot.background = element_blank(),  
    panel.background = element_blank()  
  )

#Fig2I
genes_of_interest <- c("Tbx21", "Zeb2", "Tlr7", "Myd88", "Srebf2")
expr_matrix <- GetAssayData(Bcell, slot = "data")[genes_of_interest, ]
expr_df <- as.data.frame(t(expr_matrix))
expr_df$cell <- rownames(expr_df)
meta_df <- Bcell@meta.data
meta_df$cell <- rownames(meta_df)
combined_df <- inner_join(expr_df, meta_df, by = "cell")
print(unique(combined_df$mouse))
summary_by_mouse <- combined_df %>%
  select(mouse, all_of(genes_of_interest)) %>%
  group_by(mouse) %>%
  summarise(
    across(all_of(genes_of_interest), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
summary_by_mouse$genotype <- ifelse(grepl("^FAS", summary_by_mouse$mouse), "FAS", "DKO")
summary_by_mouse <- summary_by_mouse[order(summary_by_mouse$mouse), ]
print(summary_by_mouse)
write.csv(summary_by_mouse, 
          file = "D:/FAS/Fig2/ABC_genes.csv", 
          row.names = FALSE)
