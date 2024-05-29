GetcasperHeatmap <- function(tempP, casperRData) {
  load(scDNA_RData)
  casperRData <- paste0('../005-CaSpER-UMI/',tempP,'.RData')
  load(casperRData)
  obj <- final.objects[[9]]
  casper_cnv <- obj@segments
  casper_cnv$chr <- paste0('chr',gsub('p|q','',casper_cnv$chr))
  casper_cnv[casper_cnv$states2 == 'cnloh','states2'] <- 'neut'
  casper_cnv <- casper_cnv %>% dplyr::select(ID,chr,start,end,states2)
  
  #---------step1 get signals per 1MB-----------
  bed_list <- split(casper_cnv[,2:5], casper_cnv$ID)
  casper_cnv = NULL
  for(i in 1:length(bed_list)) {
    bed = bed_list[[i]]
    gr_cnv = GRanges(seqnames = bed[, 1], ranges = IRanges(bed[, 2], bed[, 3]))
    casper_cnv = cbind(casper_cnv, average_in_window(chr_window, gr_cnv, bed[, 4],empty_v = 'neut'))
  }
  dim(casper_cnv)
  casper_cnv <- t(casper_cnv)
  rownames(casper_cnv) <- names(bed_list)
  casper_cnv[casper_cnv == 'neut'] <- '0'
  casper_cnv[casper_cnv == 'amp'] <- '1'
  casper_cnv[casper_cnv == 'del'] <- '-1'
  casper_cnv <- apply(casper_cnv, c(1, 2), function(x) ifelse(x == "0", 0, as.numeric(x)))
  casper_cnv <- apply(casper_cnv, c(1, 2), function(x) ifelse(x == "1", 1, as.numeric(x)))
  casper_cnv <- apply(casper_cnv, c(1, 2), function(x) ifelse(x == "-1", -1, as.numeric(x)))
  casper_cnv <- casper_cnv[names(cellOrder),]
  
  set.seed(123)
  hc_result <- hclust(dist(casper_cnv), method = "ward.D2")
  clusters <- cutree(hc_result, k = 2)
  C1 <- rownames(casper_cnv)[clusters == '1']
  C2 <- rownames(casper_cnv)[clusters == '2']
  
  #---------step2 define prediction------------
  casper_pre <- data.frame(Dendrogram.Group = c(rep('1',length(C1)),
                                                rep('2',length(C2))))
  rownames(casper_pre) <- c(C1,C2)
  casper_pre <- casper_pre[match(names(cellOrder),rownames(casper_pre)),,drop=F]
  
  casper_pre$scDNA <- cellOrder
  casper_pre$scRNAcasper <- 'Normal'
  hclustC1 <- casper_pre %>% filter(Dendrogram.Group == '1')
  C1_cnv <- casper_cnv[rownames(hclustC1),]
  consensusC1 <- apply(C1_cnv, 2, median)
  
  hclustC2 <- casper_pre %>% filter(Dendrogram.Group == '2')
  C2_cnv <- casper_cnv[rownames(hclustC2),]
  consensusC2 <- apply(C2_cnv, 2, median)
  
  if (sd(consensusC1) > sd(consensusC2)) {
    casper_pre[casper_pre$Dendrogram.Group == '1', 'scRNAcasper'] <- 'Tumor'
  } else {
    casper_pre[casper_pre$Dendrogram.Group == '2', 'scRNAcasper'] <- 'Tumor'
  }
  outFileName <- paste0(tempP,'_casper_byGene_summary.txt')
  sink(outFileName)
  print(casper_pre %>% count(scDNA,scRNAcasper))
  sink()
  
  outFileName <- paste0(tempP,'_casper_CNV_byGenes_byGenome.png')
  casper_rowLabel <- paste0(casper_pre$scRNAcasper,casper_pre$scDNA)
  row_ha<- rowAnnotation(casper_rowLabel = casper_rowLabel,
                         col = list(casper_rowLabel = c('NormalNormal' = '#468dc1','TumorTumor' = '#cf385e','NormalTumor'='black','TumorNormal'='black')),
                         show_legend = c('casper_rowLabel' = FALSE),
                         annotation_label = '')
  cols_anno <- c('1' = 'red', '-1' = 'blue', '0' = "white")
  
  png(outFileName,width = 800,height = 400)
  print(identical(rownames(casper_cnv),names(cellOrder)))
  ht_casper <- Heatmap(casper_cnv, name = "mat", 
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
  ht_casper <- draw(ht_casper)
  dev.off()
  outFileName <- paste0(tempP,'_casper_byGenes_ht.RDat')
  save(ht_casper,casper_cnv,casper_pre,file = outFileName)
}
