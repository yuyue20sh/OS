library(infercnv)


read_mtx <- function(data_dir, use_symbol = TRUE) {
  #
  # read matrix.mtx, barcodes.tsv, features.tsv, and meta.tsv and return a list
  #
  # Args:
  #     data_dir: str, path to the directory containing the matrix.mtx,
  #         barcodes.tsv, features.tsv, and meta.tsv files
  #
  # Returns:
  #     list, containing the matrix and meta data
  #
  matrix <- Matrix::readMM(file.path(data_dir, "matrix.mtx"))
  barcodes <- read.table(file.path(data_dir, "barcodes.tsv"), header = FALSE)
  features <- read.table(file.path(data_dir, "features.tsv"), header = FALSE)
  colnames(matrix) <- barcodes[[1]]
  if (use_symbol) {
    rownames(matrix) <- features[[2]]
  } else {
    rownames(matrix) <- features[[1]]
  }
  meta <- read.table(file.path(data_dir, "meta.tsv"), sep = "\t", header = TRUE)

  return(list(matrix, meta))
}


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
                                num_threads = 10)
  return(infercnv_obj)
}


###############################################################################

counts_dir <- "./outs/mtx/small_annotated_lvl1/"
annotations_file <-
  "/data1/hounaiqiao/yy/OS/outs/infercnv/inputs/annotation.tsv"
gene_order_file <-
  "/data1/hounaiqiao/yy/OS/outs/infercnv/inputs/hg38_gencode_v27.txt"
out_dir <- "./outs/infercnv/outputs/"

raw_counts_matrix <- as.matrix(read_mtx(counts_dir)[[1]])
ref_group_names <- c("Myeloid", "Lymphoid", "Endothelial")

infercnv_obj <- run_infercnv(raw_counts_matrix,
                             annotations_file,
                             gene_order_file,
                             ref_group_names,
                             out_dir = out_dir)

print("Done!")
