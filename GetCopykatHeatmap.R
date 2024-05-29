GetCopykatHeatmap <- function(tempP,CNVbyGeneFile,TumorNormalPredictFile) {
  copykat_cnv <- read.table(CNVbyGeneFile,header = T)
  copykat_cnv <- copykat_cnv %>% dplyr::filter(chromosome_name %in% c(1:22))
  copykat_cnv_anno <- copykat_cnv %>% dplyr::select(chromosome_name,start_position,end_position) 
  copykat_cnv_anno$chromosome_name <- paste0('chr',copykat_cnv_anno$chromosome_name)
  colnames(copykat_cnv_anno) <- c('CHR','START','END')
  copykat_cnv <- copykat_cnv[,8:ncol(copykat_cnv)]
  #--------get intersection samples and scDNA order ----------------------
  cellOrder <- cellOrder[sort(match(colnames(copykat_cnv),names(cellOrder)))]
  copykat_cnv <- copykat_cnv[,names(cellOrder)]
  #---------step1 get signals per 1MB-----------
  bed1 = cbind(copykat_cnv_anno,copykat_cnv)
  gr1 = GRanges(seqnames = bed1[, 1], ranges = IRanges(bed1[, 2], bed1[, 3]))
  copykat_cnv = average_in_window(chr_window, gr1, bed1[, -(1:3)],empty_v = 0) # change empty_v
  colnames(copykat_cnv) <- colnames(bed1[, -(1:3)])
  copykat_cnv <- t(copykat_cnv)
  #---------step2 define prediction------------
  copykat_pre <- read.table(TumorNormalPredictFile,header = T)
  copykat_pre <- copykat_pre[match(names(cellOrder),copykat_pre$cell.names),] #keep order with Cellorder
  print(identical(copykat_pre$cell.names,names(cellOrder)))
  colnames(copykat_pre)[2] <- 'scRNAcopykat'
  copykat_pre$scRNAcopykat <- if_else(grepl('diploid',copykat_pre$scRNAcopykat),'Normal','Tumor')
  copykat_pre$scDNA <- cellOrder
  outFileName <- paste0(tempP,'_copykat_byGene_summary.txt')
  sink(outFileName)
  print(copykat_pre %>% count(scDNA,scRNAcopykat))
  sink()
  
  #----------------step3 plot ht--------------
  outFileName <- paste0(tempP,'_copykat_CNV_byGenes_byGenome.png')
  copykat_rowLabel <- paste0(copykat_pre$scRNAcopykat,copykat_pre$scDNA)
  row_ha<- rowAnnotation(copykat_rowLabel = copykat_rowLabel,
                         col = list(copykat_rowLabel = c('NormalNormal' = '#468dc1','TumorTumor' = '#cf385e','NormalTumor'='black','TumorNormal'='black')),
                         show_legend = c('copykat_rowLabel' = FALSE),
                         annotation_label = '')
  cnv_vector <- as.numeric(copykat_cnv) 
  maxColor <- (summary(cnv_vector[cnv_vector>0])[6] - summary(cnv_vector[cnv_vector>0])[5])/2 + summary(cnv_vector[cnv_vector>0])[5]
  minColor <- (summary(cnv_vector[cnv_vector<0])[1] - summary(cnv_vector[cnv_vector<0])[2])/2 + summary(cnv_vector[cnv_vector<0])[2]
  
  png(outFileName,width = 800,height = 400)
  ht_copykat <-  Heatmap(copykat_cnv, name = "mat", 
                         col = colorRamp2(c(minColor, 0, maxColor), c("blue", "white", "red")), #change
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
  ht_copykat <- draw(ht_copykat)  
  dev.off()
  outFileName <- paste0(tempP,'_copykat_byGenes_ht.RDat')
  save(ht_copykat,copykat_cnv,copykat_pre,file = outFileName)
  
}





