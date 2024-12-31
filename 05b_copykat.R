# Run copykat.


library(copykat)

source("./utils/io.R")


run_copykat <- function(raw_counts_matrix, sample_name, out_dir) {
  #
  # Run copykat.
  #
  # Args:
  #     raw_counts_matrix: matrix, raw counts matrix
  #     sample_name: str, sample name
  #     out_dir: str, path to the output directory
  #
  # Returns:
  #     copykat object
  #
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


###############################################################################

data_dir <- "./outs/mtx/annotated_lvl1_split_samples/"
out_dir <- "./outs/copykat/"

dirs <- dir(data_dir)
for (d in dirs) {
  raw_counts_matrix <- as.matrix(read_mtx(file.path(data_dir, d),
                                 use_symbol = FALSE)[[1]])

  res_dir <- file.path(out_dir, d)
  if (!dir.exists(res_dir)) {
    dir.create(res_dir, recursive = TRUE)
  }

  sink(file.path(res_dir, "copykat.log"))
  run_copykat(raw_counts_matrix, d, res_dir)
  sink()
}

print("Done!")
