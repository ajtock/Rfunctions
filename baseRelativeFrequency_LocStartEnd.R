# Calculate base relative frequency around target and random loci
# Divide target loci into windows proportional to locus size
# Calculate base relative frequency for each window and
# for each nucleotide in flanking regions
# 
baseRelativeFreq <- function(targets,
                             genome,
                             mask,
                             winSize,
                             flankSize,
                             locSize,
                             outDir, 
                             targetName) {
  maskTargetOverlaps <- findOverlaps(mask, targets, ignore.strand = TRUE, select = "all")
  print(maskTargetOverlaps)
  if(length(maskTargetOverlaps) > 0) {
    targets <- targets[-subjectHits(maskTargetOverlaps)]
  }
  targets <- targets[width(targets) >= winSize]
  ranLoc <- randomizeRegions(targets, genome = genome, mask = mask, per.chromosome = TRUE, allow.overlaps = TRUE)
  
  
  a.coords <- rep(0, times = locSize)
  t.coords <- rep(0, times = locSize)
  g.coords <- rep(0, times = locSize)
  c.coords <- rep(0, times = locSize)
  ran.a.coords <- rep(0, times = locSize)
  ran.t.coords <- rep(0, times = locSize)
  ran.g.coords <- rep(0, times = locSize)
  ran.c.coords <- rep(0, times = locSize)
  target.chr.tots <- NULL
  ranLoc.chr.tots <- NULL
  for(i in 1:5) {
    print(i)        
    chr.seq <- Athaliana[[i]]
    chrTargets <- targets[seqnames(targets) == chrs[i]]
    targetStart <- start(chrTargets)-flankSize
    targetEnd <- start(chrTargets)+flankSize
    print(i)
    target.chr.tots <- c(target.chr.tots, length(chrTargets))
    print(i)
    chrRanLoc <- ranLoc[seqnames(ranLoc) == chrs[i]]
    ranLocStart <- start(chrRanLoc)-flankSize
    ranLocEnd <- start(chrRanLoc)+flankSize 
    ranLoc.chr.tots <- c(ranLoc.chr.tots, length(chrRanLoc))
    print(i)
    chr.acoords <- rep(0, times = locSize)
    chr.tcoords <- rep(0, times = locSize)
    chr.gcoords <- rep(0, times = locSize)
    chr.ccoords <- rep(0, times = locSize)
    for(j in 1:length(targetStart)) {
      acoords <- rep(0, times = locSize)
      tcoords <- rep(0, times = locSize)
      gcoords <- rep(0, times = locSize)
      ccoords <- rep(0, times = locSize)
      print(j)
      sel.seq <- unlist(strsplit(as.character(chr.seq[targetStart[j]:targetEnd[j]]), split=""))
      acoords[which(sel.seq == "A")] <- 1
      tcoords[which(sel.seq == "T")] <- 1
      gcoords[which(sel.seq == "G")] <- 1
      ccoords[which(sel.seq == "C")] <- 1 
      chr.acoords <- chr.acoords+acoords
      chr.tcoords <- chr.tcoords+tcoords
      chr.gcoords <- chr.gcoords+gcoords
      chr.ccoords <- chr.ccoords+ccoords
    }
    a.coords <- a.coords+chr.acoords
    t.coords <- t.coords+chr.tcoords
    g.coords <- g.coords+chr.gcoords
    c.coords <- c.coords+chr.ccoords
    print(i)
    chr.ran.acoords <- rep(0, times = locSize)
    chr.ran.tcoords <- rep(0, times = locSize)
    chr.ran.gcoords <- rep(0, times = locSize)
    chr.ran.ccoords <- rep(0, times = locSize)
    for(j in 1:length(ranLocStart)) {
      acoords <- rep(0, times = locSize)
      tcoords <- rep(0, times = locSize)
      gcoords <- rep(0, times = locSize)
      ccoords <- rep(0, times = locSize)
      print(j)
      sel.seq <- unlist(strsplit(as.character(chr.seq[ranLocStart[j]:ranLocEnd[j]]), split=""))
      acoords[which(sel.seq == "A")] <- 1
      tcoords[which(sel.seq == "T")] <- 1
      gcoords[which(sel.seq == "G")] <- 1
      ccoords[which(sel.seq == "C")] <- 1
      chr.ran.acoords <- chr.ran.acoords+acoords
      chr.ran.tcoords <- chr.ran.tcoords+tcoords
      chr.ran.gcoords <- chr.ran.gcoords+gcoords
      chr.ran.ccoords <- chr.ran.ccoords+ccoords
    }
    ran.a.coords <- ran.a.coords+chr.ran.acoords
    ran.t.coords <- ran.t.coords+chr.ran.tcoords
    ran.g.coords <- ran.g.coords+chr.ran.gcoords
    ran.c.coords <- ran.c.coords+chr.ran.ccoords
  }
  target.tot <- sum(target.chr.tots)
  print(target.tot)
  ranLoc.tot <- sum(ranLoc.chr.tots)
  print(ranLoc.tot)
  a.coords <- a.coords/target.tot
  t.coords <- t.coords/target.tot
  g.coords <- g.coords/target.tot
  c.coords <- c.coords/target.tot
  ran.a.coords <- ran.a.coords/ranLoc.tot
  ran.t.coords <- ran.t.coords/ranLoc.tot
  ran.g.coords <- ran.g.coords/ranLoc.tot
  ran.c.coords <- ran.c.coords/ranLoc.tot

  write.table(a.coords, file = paste0(outDir, targetName, "_targets.a.txt"))
  write.table(t.coords, file = paste0(outDir, targetName, "_targets.t.txt"))
  write.table(g.coords, file = paste0(outDir, targetName,"_targets.g.txt"))
  write.table(c.coords, file = paste0(outDir, targetName, "_targets.c.txt"))
  write.table(ran.a.coords, file = paste0(outDir, targetName, "_ranLoc.a.txt"))
  write.table(ran.t.coords, file = paste0(outDir, targetName, "_ranLoc.t.txt"))
  write.table(ran.g.coords, file = paste0(outDir, targetName, "_ranLoc.g.txt"))
  write.table(ran.c.coords, file = paste0(outDir, targetName, "_ranLoc.c.txt"))
}

