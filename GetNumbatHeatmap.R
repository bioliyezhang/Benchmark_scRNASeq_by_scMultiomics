GetNumbatHeatmap <- function(tempP, PreFileName, CNVFileName,anno.df) {
  #---------step2 define prediction------------
  numbat_pre <- read.table(PreFileName, header = T, sep = '\t') %>%
    dplyr::select(cell, clone_opt,compartment_opt)
  numbat_pre$cell <- anno.df[match(numbat_pre$cell,anno.df$ID),'CellID']
  numbat_pre <- numbat_pre[match(names(cellOrder),numbat_pre$cell),]
  print(identical(numbat_pre$cell,names(cellOrder)))
  colnames(numbat_pre)[3] <- 'scRNAnumbat'
  numbat_pre$scRNAnumbat <- if_else(grepl('normal',numbat_pre$scRNAnumbat),'Normal','Tumor')
  numbat_pre$scDNA <- cellOrder
  outFileName <- paste0(tempP,'_numbat_byGene_summary.txt')
  sink(outFileName)
  print(numbat_pre %>% count(scDNA,scRNAnumbat))
  sink()
  #---------step1 get signals per 1MB-----------
  numbat_cnv <- read.table(CNVFileName,header = T, sep = '\t') %>%
    dplyr::select(cell,CHROM,seg_start,seg_end,cnv_state_map) %>%
    distinct(cell, CHROM, seg_end, seg_start, .keep_all = T) %>%
    as.data.frame()
  numbat_cnv$CHROM <- paste0('chr',numbat_cnv$CHROM)
  
  bed_list <- split(numbat_cnv[,2:5], numbat_cnv$cell)
  numbat_cnv = NULL
  for(i in 1:length(bed_list)) {
    bed = bed_list[[i]]
    gr_cnv = GRanges(seqnames = bed[, 1], ranges = IRanges(bed[, 2], bed[, 3]))
    numbat_cnv = cbind(numbat_cnv, average_in_window(chr_window, gr_cnv, bed[, 4],empty_v = 'neu')) # change empty_v
  }
  numbat_cnv <- t(numbat_cnv)
  rownames(numbat_cnv) <- names(bed_list)
  rownames(numbat_cnv) <- anno.df[match(rownames(numbat_cnv),anno.df$ID),'CellID']
  numbat_cnv <- numbat_cnv[match(names(cellOrder),rownames(numbat_cnv)),]
  print(identical(names(cellOrder),rownames(numbat_cnv)))
  
  #----------------step3 plot ht--------------
  outFileName <- paste0(tempP,'_numbat_CNV_byGenes_byGenome.png')
  numbat_rowLabel <- paste0(numbat_pre$scRNAnumbat,numbat_pre$scDNA)
  row_ha<- rowAnnotation(numbat_rowLabel = numbat_rowLabel,
                         col = list(numbat_rowLabel = c('NormalNormal' = '#468dc1','TumorTumor' = '#cf385e','NormalTumor'='black','TumorNormal'='black')),
                         show_legend = c('numbat_rowLabel' = FALSE),
                         annotation_label = '')
  cols_anno <- c('amp' = 'darkred', 'del' = 'darkblue', 'bamp' = "salmon", 'loh' = 'darkgreen', 'bdel' = 'blue','neu'='white')
  
  png(outFileName,width = 800,height = 400)
  ht_numbat <- Heatmap(numbat_cnv, name = "mat", 
                       col = cols_anno,
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
                       use_raster = F,
                       show_heatmap_legend = F) 
  ht_numbat <- draw(ht_numbat)
  dev.off()
  outFileName <- paste0(tempP,'_numbat_byGenes_ht.RDat')
  save(ht_numbat,numbat_cnv,numbat_pre,file = outFileName)
  
}
