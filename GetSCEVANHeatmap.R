GetSCEVANHeatmap <- function(tempP,CNA_mtx,count_mtx_annot,results) {
  SCEVAN_cnv <- CNA_mtx
  SCEVAN_cnv_anno <- count_mtx_annot %>% select(seqnames,start,end)
  SCEVAN_cnv_anno$seqnames <- paste0('chr',SCEVAN_cnv_anno$seqnames)
  colnames(SCEVAN_cnv_anno) <- c('CHR','START','END')
  #--------get intersection samples and scDNA order ----------------------
  cellOrder <- cellOrder[sort(match(colnames(SCEVAN_cnv),names(cellOrder)))]
  SCEVAN_cnv <- SCEVAN_cnv[,names(cellOrder)]
  print(identical(colnames(SCEVAN_cnv),names(cellOrder)))
  #---------step1 get signals per 1MB-----------
  bed1 = cbind(SCEVAN_cnv_anno,SCEVAN_cnv)
  gr1 = GRanges(seqnames = bed1[, 1], ranges = IRanges(bed1[, 2], bed1[, 3]))
  SCEVAN_cnv = average_in_window(chr_window, gr1, bed1[, -(1:3)],empty_v = 0) # change empty_v
  colnames(SCEVAN_cnv) <- colnames(bed1[, -(1:3)])
  SCEVAN_cnv <- t(SCEVAN_cnv)
  #---------step2 define prediction------------
  SCEVAN_pre <- results
  SCEVAN_pre <- SCEVAN_pre %>% filter(class != 'filtered')
  SCEVAN_pre <- SCEVAN_pre[match(names(cellOrder),rownames(SCEVAN_pre)),]
  print(identical(rownames(SCEVAN_pre),names(cellOrder)))
  colnames(SCEVAN_pre)[1] <- 'scRNASCEVAN'
  SCEVAN_pre$scRNASCEVAN <- if_else(grepl('normal',SCEVAN_pre$scRNASCEVAN),'Normal','Tumor')
  SCEVAN_pre$scDNA <- cellOrder
  outFileName <- paste0(tempP,'_SCEVAN_byGene_summary.txt')
  sink(outFileName)
  print(SCEVAN_pre %>% count(scDNA,scRNASCEVAN))
  sink()
  #----------------step3 plot ht--------------
  outFileName <- paste0(tempP,'_SCEVAN_CNV_byGenes_byGenome.png')
  SCEVAN_rowLabel <- paste0(SCEVAN_pre$scRNASCEVAN,SCEVAN_pre$scDNA)
  row_ha<- rowAnnotation(SCEVAN_rowLabel = SCEVAN_rowLabel,
                         col = list(SCEVAN_rowLabel = c('NormalNormal' = '#468dc1','TumorTumor' = '#cf385e','NormalTumor'='black','TumorNormal'='black')),
                         show_legend = c('SCEVAN_rowLabel' = FALSE),
                         annotation_label = '')
  cnv_vector <- as.numeric(SCEVAN_cnv) 
  maxColor <- (summary(cnv_vector[cnv_vector>0])[6] - summary(cnv_vector[cnv_vector>0])[5])/2 + summary(cnv_vector[cnv_vector>0])[5]
  minColor <- (summary(cnv_vector[cnv_vector<0])[1] - summary(cnv_vector[cnv_vector<0])[2])/2 + summary(cnv_vector[cnv_vector<0])[2]
  
  png(outFileName,width = 800,height = 400)
  ht_SCEVAN <-  Heatmap(SCEVAN_cnv, name = "mat", 
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
  ht_SCEVAN <- draw(ht_SCEVAN)  
  dev.off()
  outFileName <- paste0(tempP,'_SCEVAN_byGenes_ht.RDat')
  save(ht_SCEVAN,SCEVAN_cnv,SCEVAN_pre,file = outFileName)
  
}
