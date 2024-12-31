library(infercnv)

source("./utils/io.R")


run_infercnv <- function(raw_counts_matrix,
                         annotations_file,
                         gene_order_file,
                         ref_group_names,
                         out_dir) {
  #
  # run infercnv
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
  options(scipen = 100)

  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = raw_counts_matrix,
                                       annotations_file = annotations_file,
                                       delim = "\t",
                                       gene_order_file = gene_order_file,
                                       ref_group_names = ref_group_names)

  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff = 0.1,
                                out_dir = out_dir,
                                cluster_by_groups = TRUE,
                                denoise = TRUE,
                                HMM = TRUE,
                                analysis_mode = "subclusters",
                                tumor_subcluster_partition_method =
                                  "random_trees",
                                output_format = "pdf",
                                num_threads = 20)

  return(infercnv_obj)
}


###############################################################################

data_dir <- "./outs/mtx/annotated_lvl1_split_samples/"
gene_order_file <- "./data/hg38_gencode_v27.txt"
out_dir <- "./outs/infercnv/"

dirs <- dir(data_dir)
for (d in dirs) {
  raw_counts_matrix <- as.matrix(read_mtx(file.path(data_dir, d),
                                 use_symbol = FALSE)[[1]])
  ref_group_names <- c("Monocyte/Macrophage", "Osteoclast", "Neutrophil",
                       "Mast", "T/NK", "Endothelial")
  annotations_file <- file.path(data_dir, d, "annotation.tsv")
  res_dir <- file.path(out_dir, d)
  if (!dir.exists(res_dir)) {
    dir.create(res_dir, recursive = TRUE)
  }
  sink(file.path(res_dir, "infercnv.log"))
  infercnv_obj <- run_infercnv(raw_counts_matrix,
                               annotations_file,
                               gene_order_file,
                               ref_group_names,
                               out_dir = res_dir)
  sink()
}

print("Done!")
