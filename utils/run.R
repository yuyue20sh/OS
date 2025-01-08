run_infercnv <- function(raw_counts_matrix,
                         annotations_file,
                         gene_order_file,
                         ref_group_names,
                         out_dir) {
  #
  # Run infercnv.
  #
  # Args:
  #     raw_counts_matrix: matrix, raw counts matrix
  #     annotations_file: str, path to the annotations file
  #     gene_order_file: str, path to the gene order file
  #     ref_group_names: list, reference group names
  #     out_dir: str, path to the output directory
  #
  # Returns:
  #     infercnv object
  #
  library(infercnv)
  options(scipen = 100)

  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = raw_counts_matrix,
                                       annotations_file = annotations_file,
                                       delim = "\t",
                                       gene_order_file = gene_order_file,
                                       ref_group_names = ref_group_names)

  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff = 0.1,
                                out_dir = out_dir,
                                cluster_by_groups = FALSE,
                                denoise = TRUE,
                                HMM = TRUE,
                                analysis_mode = "subclusters",
                                tumor_subcluster_partition_method =
                                  "random_trees",
                                output_format = "pdf",
                                num_threads = 1)

  return(infercnv_obj)
}


run_copykat <- function(raw_counts_matrix, sample_name, out_dir) {
  #
  # Run copykat.
  #
  # Args:
  #     raw_counts_matrix: matrix, raw counts matrix
  #     sample_name: str, sample name
  #     out_dir: str, path to the output directory
  #
  library(copykat)
  copykat_res <- copykat(rawmat = raw_counts_matrix,
                         id.type = "S",
                         ngene.chr = 5,
                         LOW.DR = 0.05,
                         UP.DR = 0.2,
                         win.size = 25,
                         KS.cut = 0.1,
                         distance = "euclidean",
                         sam.name = paste0(out_dir, "/", sample_name),
                         norm.cell.names = "",
                         plot.genes = "TRUE",
                         genome = "hg20",
                         n.cores = 20)

  # save
  saveRDS(copykat_res, file = file.path(out_dir,
                                        paste0(sample_name, ".copykat.rds")))

  return(copykat_res)
}
