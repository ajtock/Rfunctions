# Function to create coverage matrices for target loci and random loci (incl. flanking regions)
## and to calculate mean levels per window
covMatrix <- function(signal, target, ranLoc, targetSize, flankSize, winSize, outDFCM, x) {
  # target loci
  set.seed(2840)
  mat1_smoothed <- normalizeToMatrix(signal, target, value_column = "coverage",
                                     extend = flankSize, mean_mode = "absolute", w = winSize,
                                     background = 0, smooth = TRUE,
                                     include_target = TRUE, target_ratio = targetSize/(targetSize+(flankSize*2)))
  print(mat1_smoothed)
  print(length(mat1_smoothed))
  mat1_smoothed_DF <- data.frame(mat1_smoothed)
  mat1_smoothed_DF_colMeans <- as.vector(colMeans(mat1_smoothed_DF))
  write.table(mat1_smoothed_DF_colMeans, file = outDFCM[[x]][[1]])

  # random loci
  set.seed(8472)
  mat2_smoothed <- normalizeToMatrix(signal, ranLoc, value_column = "coverage",
                                     extend = flankSize, mean_mode = "absolute", w = winSize,
                                     background = 0, smooth = TRUE,
                                     include_target = TRUE, target_ratio = targetSize/(targetSize+(flankSize*2)))
  print(mat2_smoothed)
  print(length(mat2_smoothed))
  mat2_smoothed_DF <- data.frame(mat2_smoothed)
  mat2_smoothed_DF_colMeans <- as.vector(colMeans(mat2_smoothed_DF))
  write.table(mat2_smoothed_DF_colMeans, file = outDFCM[[x]][[2]])
}

# Function to create DNA methylation matrices for target loci and random loci (incl. flanking regions)
## and to calculate mean levels per window
DNAmethMatrix <- function(signal, target, ranLoc, targetSize, flankSize, winSize, DNAmethOutDFCM, x) {
  # target loci
  set.seed(2840)
  methMat1_smoothed <- normalizeToMatrix(signal, target, value_column = "coverage",
                                         extend = flankSize, mean_mode = "absolute", w = winSize,
                                         background = NA, smooth = TRUE,
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

  write.table(methMat1_smoothed_DF_colMeans, file = DNAmethOutDFCM[[x]][[1]])

  # random loci
  set.seed(8472)
  methMat2_smoothed <- normalizeToMatrix(signal, ranLoc, value_column = "coverage",
                                         extend = flankSize, mean_mode = "absolute", w = winSize,
                                         background = NA, smooth = TRUE,
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

  write.table(methMat2_smoothed_DF_colMeans, file = DNAmethOutDFCM[[x]][[2]])
}

