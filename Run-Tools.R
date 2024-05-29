
#----------------inferCNV----------------
infercnv_obj = CreateInfercnvObject(delim = '\t',
                                      raw_counts_matrix = epi_RNA,
                                      annotations_file = tempFile,
                                      gene_order_file = 'hg38_gencode_v27.txt',
                                      ref_group_names = NULL)
tempFile <- paste0(tempP,'_Epi_HMMFalse')
infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, 
                               k_obs_groups = 2,
                               cluster_by_groups = F,
                               out_dir=tempFile,
                               denoise=TRUE,
                               HMM=F,
                               tumor_subcluster_partition_method='random_trees',
                               tumor_subcluster_pval=0.1,
                               num_threads = 1,
                               output_format = 'pdf')

#----------------CopyKAT----------------
copykat.test <- copykat(rawmat=temp_RNA, 
                          id.type="S", 
                          cell.line = "no",
                          ngene.chr=1, 
                          LOW.DR = 0.05,
                          UP.DR = 0.1,
                          win.size=25, 
                          norm.cell.names = "",
                          KS.cut=0.1, 
                          sam.name=tempP, 
                          distance="euclidean",
                          output.seg="FLASE", 
                          plot.genes="TRUE", 
                          genome="hg20",
                          n.cores=1)

#----------------SCEVAN----------------
results <- pipelineCNA(temp_RNA, 
                         sample = tempP,
                         par_cores = 1,
                         norm_cell = NULL,
                         SUBCLONES = F,
                         beta_vega = 0.5,
                         ClonalCN = F,
                         plotTree = F,
                         AdditionalGeneSets = NULL,
                         SCEVANsignatures = TRUE,
                         organism = "human")

#----------------Numbat----------------
out = run_numbat( count_mat = temp_RNA,
                    lambdas_ref = temp_ref, #
                    df_allele = temp_allel,
                    genome = "hg38",
                    t = 1e-5,
                    ncores = 5,
                    plot = TRUE,
                    gamma = 5,
                    min_genes = 0,min_cells = 10,
                    out_dir = tempP)

#----------------CaSpER----------------
object <- CreateCasperObject(raw.data=temp_RNA,
                               annotation=annotation,
                               control.sample.ids = NormCells,
                               cytoband=cytoband_hg38,
                               loh.name.mapping=temp.loh.name.mapping, 
                               cnv.scale=3, 
                               loh.scale=3,
                               method="iterative",
                               loh=temp.loh, 
                               sequencing.type="single-cell",
                               log.transformed = F,
                               genomeVersion = 'hg38')
  final.objects <- runCaSpER(object,
                             removeCentromere=T, 
                             cytoband=cytoband_hg38, 
                             method="iterative")