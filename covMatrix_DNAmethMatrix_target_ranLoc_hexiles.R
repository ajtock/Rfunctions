# Function to create coverage matrices for target loci and random loci (incl. flanking regions)
## and to calculate mean levels per window
covMatrixHexile <- function(signal, targetsHex, targetSize, flankSize, winSize, outDFCM, x, y) {
  # target loci
  target <- read.table(targetsHex[y], header = T)
  target <- target[,1:5]
  colnames(target) <- c("chr", "start", "end", "strand", "gene_model")
  targetGR <- GRanges(seqnames = target$chr, ranges = IRanges(start = target$start, end = target$end),
                      strand = target$strand, gene_model = target$gene_model)

  set.seed(2840)
  mat1_smoothed <- normalizeToMatrix(signal, targetGR, value_column = "coverage",
                                     extend = flankSize, mean_mode = "absolute", w = winSize,
                                     empty_value = 0, smooth = TRUE,
                                     include_target = TRUE, target_ratio = targetSize/(targetSize+(flankSize*2)))
  print(mat1_smoothed)
  print(length(mat1_smoothed))
  mat1_smoothed_DF <- data.frame(mat1_smoothed)
  mat1_smoothed_DF_colMeans <- as.vector(colMeans(mat1_smoothed_DF))
  write.table(mat1_smoothed_DF_colMeans, file = outDFCM[[x]][[1]][y])

  # random loci
  # Generate GRanges object containing random loci of same number and size distribution as targetGR
  ranLocGR <- randomizeRegions(targetGR, genome = genome, mask = mask, per.chromosome = TRUE, allow.overlaps = TRUE)

  set.seed(8472)
  mat2_smoothed <- normalizeToMatrix(signal, ranLocGR, value_column = "coverage",
                                     extend = flankSize, mean_mode = "absolute", w = winSize,
                                     empty_value = 0, smooth = TRUE,
                                     include_target = TRUE, target_ratio = targetSize/(targetSize+(flankSize*2)))
  print(mat2_smoothed)
  print(length(mat2_smoothed))
  mat2_smoothed_DF <- data.frame(mat2_smoothed)
  mat2_smoothed_DF_colMeans <- as.vector(colMeans(mat2_smoothed_DF))
  write.table(mat2_smoothed_DF_colMeans, file = outDFCM[[x]][[2]][y])
}

# Function to create DNA methylation matrices for target loci and random loci (incl. flanking regions)
## and to calculate mean levels per window
DNAmethMatrixHexile <- function(signal, targetsHex, targetSize, flankSize, winSize, DNAmethOutDFCM, x, y) {
  # target loci
  target <- read.table(targetsHex[y], header = T)
  target <- target[,1:5]
  colnames(target) <- c("chr", "start", "end", "strand", "gene_model")
  targetGR <- GRanges(seqnames = target$chr, ranges = IRanges(start = target$start, end = target$end),
                      strand = target$strand, gene_model = target$gene_model)

  set.seed(2840)
  methMat1_smoothed <- normalizeToMatrix(signal, targetGR, value_column = "coverage",
                                         extend = flankSize, mean_mode = "absolute", w = winSize,
                                         empty_value = NA, smooth = TRUE,
                                         include_target = TRUE, target_ratio = targetSize/(targetSize+(flankSize*2)))
  print(methMat1_smoothed)
  print(length(methMat1_smoothed))
  print("methMat1_smoothed rows = ")
  print(length(methMat1_smoothed)/round((targetSize/winSize)+((flankSize*2)/winSize)))
  methMat1_smoothed_failed_rows <- attr(methMat1_smoothed, "failed_rows")
  print("methMat1_smoothed failed rows = ")
  print(length(methMat1_smoothed_failed_rows))
  if(is.null(methMat1_smoothed_failed_rows) == FALSE) {
    methMat1_smoothed <- methMat1_smoothed[-methMat1_smoothed_failed_rows,]
  }
  print(methMat1_smoothed)
  print(length(methMat1_smoothed))
  print("methMat1_smoothed rows less failed rows = ")
  print(length(methMat1_smoothed)/round((targetSize/winSize)+((flankSize*2)/winSize)))
  methMat1_smoothed_DF <- data.frame(methMat1_smoothed)
  methMat1_smoothed_DF_colMeans <- as.vector(colMeans(methMat1_smoothed_DF))

  write.table(methMat1_smoothed_DF_colMeans, file = DNAmethOutDFCM[[x]][[1]][y])

  # random loci
  # Generate GRanges object containing random loci of same number and size distribution as targetGR
  ranLocGR <- randomizeRegions(targetGR, genome = genome, mask = mask, per.chromosome = TRUE, allow.overlaps = TRUE)

  set.seed(8472)
  methMat2_smoothed <- normalizeToMatrix(signal, ranLocGR, value_column = "coverage",
                                         extend = flankSize, mean_mode = "absolute", w = winSize,
                                         empty_value = NA, smooth = TRUE,
                                         include_target = TRUE, target_ratio = targetSize/(targetSize+(flankSize*2)))
  print(methMat2_smoothed)
  print(length(methMat2_smoothed))
  print("methMat2_smoothed rows = ")
  print(length(methMat2_smoothed)/round((targetSize/winSize)+((flankSize*2)/winSize)))
  methMat2_smoothed_failed_rows <- attr(methMat2_smoothed, "failed_rows")
  print("methMat2_smoothed failed rows = ")
  print(length(methMat2_smoothed_failed_rows))
  if(is.null(methMat2_smoothed_failed_rows) == FALSE) {
    methMat2_smoothed <- methMat2_smoothed[-methMat2_smoothed_failed_rows,]
  }
  print(methMat2_smoothed)
  print(length(methMat2_smoothed))
  print("methMat2_smoothed rows less failed rows = ")
  print(length(methMat2_smoothed)/round((targetSize/winSize)+((flankSize*2)/winSize)))
  methMat2_smoothed_DF <- data.frame(methMat2_smoothed)
  methMat2_smoothed_DF_colMeans <- as.vector(colMeans(methMat2_smoothed_DF))

  write.table(methMat2_smoothed_DF_colMeans, file = DNAmethOutDFCM[[x]][[2]][y])
}

