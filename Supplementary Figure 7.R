#Supplementary Figure 7
library(ggplot2)
library(Seurat)

#FigS7A
#Fig1B
feature_genes<-c("MS4A1","CD27","TCL1A","CD83","ITGAX","TBX21")
output_dir <- "D:/FAS/Fig7/"
dir.create(output_dir, showWarnings = FALSE)
for (gene in feature_genes) {
  if (!gene %in% rownames(SLE)) {
    warning(paste("Gene", gene, "not found. Skipping..."))
    next
  }
  
  p <- FeaturePlot(SLE, features = gene, order = TRUE, pt.size = 1.2) +
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


#FigS7B
feature_genes<-c("MS4A1","CD27","TCL1A","TNFRSF13B","ITGAX","TBX21")
output_dir <- "D:/FAS/Fig7/"
dir.create(output_dir, showWarnings = FALSE)
for (gene in feature_genes) {
  if (!gene %in% rownames(pSS)) {
    warning(paste("Gene", gene, "not found. Skipping..."))
    next
  }
  
  p <- FeaturePlot(pSS, features = gene, order = TRUE, pt.size = 1.2) +
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



#FigS7C
feature_genes<-c("MS4A1","CD27","TCL1A","TNFRSF13B","ITGAX","TBX21")
output_dir <- "D:/FAS/Fig7/"
dir.create(output_dir, showWarnings = FALSE)
for (gene in feature_genes) {
  if (!gene %in% rownames(MS)) {
    warning(paste("Gene", gene, "not found. Skipping..."))
    next
  }
  
  p <- FeaturePlot(MS, features = gene, order = TRUE, pt.size = 1.2) +
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
