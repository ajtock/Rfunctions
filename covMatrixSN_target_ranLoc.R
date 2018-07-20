# Function to create coverage matrices for target loci and random loci
# (single-nucleotide positions) and flanking regions and to calculate mean levels per window
covMatrixSN <- function(signal, target, ranLoc, flankSize, winSize, outDFCM, x) {
  # target loci
  set.seed(2840)
  mat1_smoothed <- normalizeToMatrix(signal, target, value_column = "coverage",
                                     extend = flankSize, mean_mode = "absolute", w = winSize,
                                     empty_value = 0, smooth = TRUE,
                                     include_target = FALSE)
  print(mat1_smoothed)
  print(length(mat1_smoothed))
  mat1_smoothed_DF <- data.frame(mat1_smoothed)
  mat1_smoothed_DF_colMeans <- as.vector(colMeans(mat1_smoothed_DF))
  write.table(mat1_smoothed_DF_colMeans, file = outDFCM[[x]][[1]])

  # random loci
  set.seed(8472)
  mat2_smoothed <- normalizeToMatrix(signal, ranLoc, value_column = "coverage",
                                     extend = flankSize, mean_mode = "absolute", w = winSize,
                                     empty_value = 0, smooth = TRUE,
                                     include_target = FALSE)
  print(mat2_smoothed)
  print(length(mat2_smoothed))
  mat2_smoothed_DF <- data.frame(mat2_smoothed)
  mat2_smoothed_DF_colMeans <- as.vector(colMeans(mat2_smoothed_DF))
  write.table(mat2_smoothed_DF_colMeans, file = outDFCM[[x]][[2]])
}

