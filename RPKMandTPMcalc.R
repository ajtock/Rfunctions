# Function to calculate RPKM and TPM,
# and generate feature GRanges object with these values appended
RPKMandTPMcalcSE <- function(bamFile, libName, features, featureName, outDir) {
  # Load BAM and create RangedData object
  lib_readGAlignments <- readGAlignments(bamFile)
  lib_ranges <- ranges(lib_readGAlignments)
  lib_chrs <- as.data.frame(seqnames(lib_readGAlignments))
  lib_chrs <- as.character(lib_chrs[,1])
  lib_strand <- as.character(strand(lib_readGAlignments))
  lib_ranged <- RangedData(ranges = lib_ranges, space = lib_chrs, strand = lib_strand)
  save(lib_ranged, file = paste0(outDir, libName, "_RangedData.RData"))

  # Calculate library size
  chr_size <- NULL
  for(j in 1:5) {
    print(j)
    chr_lib_ranged <- lib_ranged[j]
    chr_size <- c(chr_size, length(space(chr_lib_ranged)))
  }
  lib_size <- sum(chr_size)
  # Calculate "per million" scaling factor
  RPM_scaling_factor <- lib_size/1e+06

  # Calculate RPM and RPKM for each feature
  lib_rangedGR <- as(lib_ranged, "GRanges")

  feature_reads <- countOverlaps(features, lib_rangedGR)
  feature_RPM <- feature_reads/RPM_scaling_factor
  feature_RPKM <- feature_RPM/(width(features)/1e+03)

  # Calculate TPM for each feature
  feature_RPK <- feature_reads/(width(features)/1e+03)
  RPKPM_scaling_factor <- sum(feature_RPK)/1e+06
  feature_TPM <- feature_RPK/RPKPM_scaling_factor

  # Append feature_TPM and feature_RPKM as extra fields in features GRanges object
  featuresGR <- GRanges(seqnames = seqnames(features), ranges = ranges(features), strand = strand(features), feature_model = features$gene_model, feature_TPM, feature_RPKM)
  save(featuresGR, file = paste0(outDir, libName, "_TPM_RPKM_at_", featureName, ".RData"))
  featuresGRsortedTPM <- sort(featuresGR, by = ~ feature_TPM, decreasing = T)
  save(featuresGRsortedTPM, file = paste0(outDir, libName, "_TPM_RPKM_at_", featureName, "_sorted_by_decreasing_TPM.RData"))
}

# Function to calculate RPKM and TPM,
# and generate feature GRanges object with these values appended
RPKMandTPMcalcPE <- function(bamFile, libName, features, featureName, outDir) {
  # Load BAM and create RangedData object
  lib_readGAlignmentPairs <- readGAlignmentPairs(bamFile)
  lib_ranges <- ranges(lib_readGAlignmentPairs)
  lib_chrs <- as.data.frame(seqnames(lib_readGAlignmentPairs))
  lib_chrs <- as.character(lib_chrs[,1])
  lib_strand <- as.character(strand(lib_readGAlignmentPairs))
  lib_ranged <- RangedData(ranges = lib_ranges, space = lib_chrs, strand = lib_strand)
  save(lib_ranged, file = paste0(outDir, libName, "_RangedData.RData"))

  # Calculate library size
  chr_size <- NULL
  for(j in 1:5) {
    print(j)
    chr_lib_ranged <- lib_ranged[j]
    chr_size <- c(chr_size, length(space(chr_lib_ranged)))
  }
  lib_size <- sum(chr_size)
  # Calculate "per million" scaling factor
  RPM_scaling_factor <- lib_size/1e+06

  # Calculate RPM and RPKM for each feature
  lib_rangedGR <- as(lib_ranged, "GRanges")

  feature_reads <- countOverlaps(features, lib_rangedGR)
  feature_RPM <- feature_reads/RPM_scaling_factor
  feature_RPKM <- feature_RPM/(width(features)/1e+03)

  # Calculate TPM for each feature
  feature_RPK <- feature_reads/(width(features)/1e+03)
  RPKPM_scaling_factor <- sum(feature_RPK)/1e+06
  feature_TPM <- feature_RPK/RPKPM_scaling_factor

  # Append feature_TPM and feature_RPKM as extra fields in features GRanges object
  featuresGR <- GRanges(seqnames = seqnames(features), ranges = ranges(features), strand = strand(features), feature_model = features$gene_model, feature_TPM, feature_RPKM)
  save(featuresGR, file = paste0(outDir, libName, "_TPM_RPKM_at_", featureName, ".RData"))
  featuresGRsortedTPM <- sort(featuresGR, by = ~ feature_TPM, decreasing = T)
  save(featuresGRsortedTPM, file = paste0(outDir, libName, "_TPM_RPKM_at_", featureName, "_sorted_by_decreasing_TPM.RData"))
}

