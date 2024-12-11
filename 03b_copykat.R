library(copykat)


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


###############################################################################

data_dir <- "./outs/mtx/concatenated/"

for (d in list.dirs(data_dir, recursive = FALSE)) {
  # read data
  id <- basename(d)
  print(id)
  data <- read_mtx(d)
  matrix <- data[[1]]
  meta <- data[[2]]
  rownames(matrix) <- make.unique(rownames(matrix))

  # run copykat
  sam_dir <- file.path("./outs/copykat/", id)
  if (!dir.exists(sam_dir)) {
    dir.create(sam_dir, recursive = TRUE)
  }
  copykat_res <-
    copykat(rawmat = as.matrix(matrix),
            id.type = "S",
            ngene.chr = 5,
            LOW.DR = 0.01,
            win.size = 25,
            KS.cut = 0.1,
            distance = "euclidean",
            sam.name = paste0(sam_dir, "/", id),
            norm.cell.names = "",
            plot.genes = "TRUE",
            genome = "hg20",
            n.cores = 20)

  # save
  saveRDS(copykat_res, file = paste0(sam_dir, "/", id, ".copykat.rds"))
}

print("Done!")

# "low confidence in classification" (as well as 'WARNING! NOT CONVERGENT!'):
# BC3, DLJM44, WBXM16, WBXM16_2, WYQM12 and XZHM13
# "WARNING! NOT CONVERGENT!": ZCLM12
