GetinferCNVHeatmap <- function(tempP,infercnvRData,PreFileName) {
  inferCNV_obj <- readRDS(infercnvRData)
  inferCNV_cnv <- inferCNV_obj@expr.data
  gene_order <- inferCNV_obj@gene_order
  colnames(gene_order) <- c('CHR','START','END')
  inferCNV_cnv <- inferCNV_cnv[,match(names(cellOrder),colnames(inferCNV_cnv))] # ordered by chromosome and Normal+Tumor
  print(identical(names(cellOrder),colnames(inferCNV_cnv)))
  
  #---------step1 get signals per 1MB-----------
  bed1 = cbind(gene_order,inferCNV_cnv)
  gr1 = GRanges(seqnames = bed1[, 1], ranges = IRanges(bed1[, 2], bed1[, 3]))
  inferCNV_cnv = average_in_window(chr_window, gr1, bed1[, -(1:3)],empty_v = 1) # change empty_v
  colnames(inferCNV_cnv) <- colnames(bed1[, -(1:3)])
  inferCNV_cnv <- t(inferCNV_cnv)
  
  #---------step2 define prediction------------
  inferCNV_pre <- read.table(PreFileName)
  inferCNV_pre <- inferCNV_pre[match(names(cellOrder),rownames(inferCNV_pre)),]
  inferCNV_pre$scDNA <- cellOrder
  inferCNV_pre$scRNAinferCNV <- 'Normal'
  
  hclustC1 <- inferCNV_pre %>% filter(Dendrogram.Group == '1')
  C1_cnv <- inferCNV_cnv[rownames(hclustC1),]
  consensusC1 <- apply(C1_cnv, 2, median)
  
  hclustC2 <- inferCNV_pre %>% filter(Dendrogram.Group == '2')
  C2_cnv <- inferCNV_cnv[rownames(hclustC2),]
  consensusC2 <- apply(C2_cnv, 2, median)
  
  if (sd(consensusC1) > sd(consensusC2)) {
    inferCNV_pre[inferCNV_pre$Dendrogram.Group == '1', 'scRNAinferCNV'] <- 'Tumor'
  } else {
    inferCNV_pre[inferCNV_pre$Dendrogram.Group == '2', 'scRNAinferCNV'] <- 'Tumor'
  }
  
  outFileName <- paste0(tempP,'_inferCNV_byGene_summary.txt')
  sink(outFileName)
  print(inferCNV_pre %>% count(scDNA,scRNAinferCNV))
  sink()
  #----------------step3 plot ht--------------
  outFileName <- paste0(tempP,'_inferCNV_CNV_byGenes_byGenome.png')
  inferCNV_rowLabel <- paste0(inferCNV_pre$scRNAinferCNV,inferCNV_pre$scDNA)
  row_ha<- rowAnnotation(inferCNV_rowLabel = inferCNV_rowLabel,
                         col = list(inferCNV_rowLabel = c('NormalNormal' = '#468dc1','TumorTumor' = '#cf385e','NormalTumor'='black','TumorNormal'='black')),
                         show_legend = c('inferCNV_rowLabel' = FALSE),
                         annotation_label = '')
  cnv_vector <- as.numeric(inferCNV_cnv) 
  maxColor <- (summary(cnv_vector[cnv_vector>1])[6] - summary(cnv_vector[cnv_vector>1])[5])/2 + summary(cnv_vector[cnv_vector>1])[5]
  minColor <- (summary(cnv_vector[cnv_vector<1])[1] - summary(cnv_vector[cnv_vector<1])[2])/2 + summary(cnv_vector[cnv_vector<1])[2]
  
  png(outFileName,width = 800,height = 400)
  ht_inferCNV <-  Heatmap(inferCNV_cnv, name = "mat", 
                          col = colorRamp2(c(minColor, 1, maxColor), c("blue", "white", "red")),
                          left_annotation = row_ha, 
                          row_title = NULL,
                          row_split = cellOrder, 
                          cluster_rows = F, 
                          show_row_names = F,
                          show_row_dend = F,
                          row_dend_width = unit(0.5, "cm"),
                          cluster_columns = F,
                          column_split = chr,
                          column_title_gp = gpar(fontsize = 8),
                          border = TRUE,
                          use_raster = T,
                          show_heatmap_legend = F) 
  ht_inferCNV <- draw(ht_inferCNV)  
  dev.off()
  outFileName <- paste0(tempP,'_inferCNV_byGenes_ht.RDat')
  save(ht_inferCNV,inferCNV_cnv,inferCNV_pre,file = outFileName)
  
}
