GetscDNAHeatmap <- function(tempP) {
  temp_cnv <- cnv_mat[grep(tempP,rownames(cnv_mat)),]
  temp_meta <- epi_meta[match(rownames(temp_cnv),epi_meta$CB),]
  identical(temp_meta$CB,rownames(temp_cnv))
  celltype <- temp_meta$CellType
  row_ha<- rowAnnotation(celltype = celltype, 
                         col = list(celltype = c('Normal' = '#468dc1','Tumor' = '#cf385e')),
                         show_legend = c('celltype' = FALSE),
                         annotation_label = "")
  
  outFileName <- paste0('scDNA_CNV_',tempP,'_byGenome_dend.png')
  cnv_vector <- as.numeric(temp_cnv) 
  maxColor <- (summary(cnv_vector[cnv_vector>2])[6] - summary(cnv_vector[cnv_vector>2])[5])/2 + summary(cnv_vector[cnv_vector>2])[5]
  minColor <- (summary(cnv_vector[cnv_vector<2])[1] - summary(cnv_vector[cnv_vector<2])[2])/2 + summary(cnv_vector[cnv_vector<2])[2]
  png(outFileName,width = 800,height = 400)
  ht <- Heatmap(temp_cnv, name = "mat", 
                col = colorRamp2(c(minColor, 2, maxColor), c("blue", "white", "red")),
                left_annotation = row_ha, 
                row_title = NULL,
                row_split = celltype, 
                cluster_rows = T, 
                show_row_names = F,
                show_row_dend = F,
                row_dend_width = unit(0.5, "cm"),
                cluster_columns = F,
                column_split = chr,
                column_title_gp = gpar(fontsize = 8),
                border = TRUE,
                use_raster = T,
                show_heatmap_legend = F) 
  
  ht <- draw(ht)  
  dev.off()
  
  NormalCells <- rownames(temp_cnv)[row_order(ht)$Normal]
  TumorCells <- rownames(temp_cnv)[row_order(ht)$Tumor]
  cellOrder <- c(rep('Normal',length(NormalCells)),rep('Tumor',length(TumorCells)))
  names(cellOrder) <- c(NormalCells,TumorCells)
  outFileName <- paste0('scDNA_CNV_',tempP,'_ht.RData')
  save(temp_cnv,temp_meta,row_ha,ht,cellOrder,file = outFileName)
}