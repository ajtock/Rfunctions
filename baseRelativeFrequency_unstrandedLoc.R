## Define "rounding in commerce" function so that numbers ending .5 round up (not
## to the nearest even digit as in the R round() function)
## See http://andland.github.io/blog/2012/06/15/rounding-in-r/
#cround <- function(x,n) {
#  vorz = sign(x)
#  z = abs(x)*10^n
#  z = z + 0.5
#  z = trunc(z)
#  z = z/10^n
#  z*vorz
#}

# Calculate base relative frequency around taget and random loci
baseRelativeFreqUnstranded <- function(targets, genome, mask, flankSize, locSize, outDir, targetName) {
  maskTargetOverlaps <- findOverlaps(mask, targets, ignore.strand = TRUE, select = "all")
  print(maskTargetOverlaps)
  if(length(maskTargetOverlaps) > 0) {
    targets <- targets[-subjectHits(maskTargetOverlaps)]
  }
  set.seed(374592)
  ranLoc <- randomizeRegions(targets, genome = genome, mask = mask, per.chromosome = TRUE, allow.overlaps = TRUE)
  a.coords.targetPlus <- rep(0, times = locSize)
  t.coords.targetPlus <- rep(0, times = locSize)
  g.coords.targetPlus <- rep(0, times = locSize)
  c.coords.targetPlus <- rep(0, times = locSize)
  a.coords.ranLocPlus <- rep(0, times = locSize)
  t.coords.ranLocPlus <- rep(0, times = locSize)
  g.coords.ranLocPlus <- rep(0, times = locSize)
  c.coords.ranLocPlus <- rep(0, times = locSize)

  a.coords.targetMinus <- rep(0, times = locSize)
  t.coords.targetMinus <- rep(0, times = locSize)
  g.coords.targetMinus <- rep(0, times = locSize)
  c.coords.targetMinus <- rep(0, times = locSize)
  a.coords.ranLocMinus <- rep(0, times = locSize)
  t.coords.ranLocMinus <- rep(0, times = locSize)
  g.coords.ranLocMinus <- rep(0, times = locSize)
  c.coords.ranLocMinus <- rep(0, times = locSize)

  targetPlus.chr.tots <- NULL
  ranLocPlus.chr.tots <- NULL
  targetMinus.chr.tots <- NULL
  ranLocMinus.chr.tots <- NULL
  for(i in 1:5) {
    chr.seq <- Athaliana[[i]]
    chrTarget <- targets[seqnames(targets) == chrs[i]]
    halfTargetLength <- round(length(chrTarget)/2)
    chrTargetPlus <- chrTarget[1:halfTargetLength]
    chrTargetMinus <- chrTarget[(halfTargetLength+1):length(chrTarget)]
    targetPlus.chr.tots <- c(targetPlus.chr.tots, length(chrTargetPlus))
    targetMinus.chr.tots <- c(targetMinus.chr.tots, length(chrTargetMinus))
    targetPlusStart <- start(chrTargetPlus)-flankSize
    targetPlusEnd <- start(chrTargetPlus)+flankSize
    targetMinusStart <- start(chrTargetMinus)-flankSize
    targetMinusEnd <- start(chrTargetMinus)+flankSize
    print(i)

    chrRanLoc <- ranLoc[seqnames(ranLoc) == chrs[i]]
    halfRanLocLength <- round(length(chrRanLoc)/2)
    chrRanLocPlus <- chrRanLoc[1:halfRanLocLength]
    chrRanLocMinus <- chrRanLoc[(halfRanLocLength+1):length(chrRanLoc)]
    ranLocPlus.chr.tots <- c(ranLocPlus.chr.tots, length(chrRanLocPlus))
    ranLocMinus.chr.tots <- c(ranLocMinus.chr.tots, length(chrRanLocMinus))
    ranLocPlusStart <- start(chrRanLocPlus)-flankSize
    ranLocPlusEnd <- start(chrRanLocPlus)+flankSize
    ranLocMinusStart <- start(chrRanLocMinus)-flankSize
    ranLocMinusEnd <- start(chrRanLocMinus)+flankSize
    print(i)

    # "Plus" strand loci
    chr.acoords.targetPlus <- rep(0, times = locSize)
    chr.tcoords.targetPlus <- rep(0, times = locSize)
    chr.gcoords.targetPlus <- rep(0, times = locSize)
    chr.ccoords.targetPlus <- rep(0, times = locSize)
    chr.acoords.ranLocPlus <- rep(0, times = locSize)
    chr.tcoords.ranLocPlus <- rep(0, times = locSize)
    chr.gcoords.ranLocPlus <- rep(0, times = locSize)
    chr.ccoords.ranLocPlus <- rep(0, times = locSize)
    print(i)
    for(j in 1:length(targetPlusStart)) {
      print(i)
      acoords.targetPlus <- rep(0, times = locSize)
      tcoords.targetPlus <- rep(0, times = locSize)
      gcoords.targetPlus <- rep(0, times = locSize)
      ccoords.targetPlus <- rep(0, times = locSize)
      print(j)

      sel.seq.targetPlus <- unlist(strsplit(as.character(chr.seq[targetPlusStart[j]:targetPlusEnd[j]]), split=""))
      acoords.targetPlus[which(sel.seq.targetPlus == "A")] <- 1
      tcoords.targetPlus[which(sel.seq.targetPlus == "T")] <- 1
      gcoords.targetPlus[which(sel.seq.targetPlus == "G")] <- 1
      ccoords.targetPlus[which(sel.seq.targetPlus == "C")] <- 1 
      chr.acoords.targetPlus <- chr.acoords.targetPlus+acoords.targetPlus
      chr.tcoords.targetPlus <- chr.tcoords.targetPlus+tcoords.targetPlus
      chr.gcoords.targetPlus <- chr.gcoords.targetPlus+gcoords.targetPlus
      chr.ccoords.targetPlus <- chr.ccoords.targetPlus+ccoords.targetPlus
      print(j)

      acoords.ranLocPlus <- rep(0, times = locSize)
      tcoords.ranLocPlus <- rep(0, times = locSize)
      gcoords.ranLocPlus <- rep(0, times = locSize)
      ccoords.ranLocPlus <- rep(0, times = locSize)
      print(j)

      sel.seq.ranLocPlus <- unlist(strsplit(as.character(chr.seq[ranLocPlusStart[j]:ranLocPlusEnd[j]]), split=""))
      acoords.ranLocPlus[which(sel.seq.ranLocPlus == "A")] <- 1
      tcoords.ranLocPlus[which(sel.seq.ranLocPlus == "T")] <- 1
      gcoords.ranLocPlus[which(sel.seq.ranLocPlus == "G")] <- 1
      ccoords.ranLocPlus[which(sel.seq.ranLocPlus == "C")] <- 1 
      chr.acoords.ranLocPlus <- chr.acoords.ranLocPlus+acoords.ranLocPlus
      chr.tcoords.ranLocPlus <- chr.tcoords.ranLocPlus+tcoords.ranLocPlus
      chr.gcoords.ranLocPlus <- chr.gcoords.ranLocPlus+gcoords.ranLocPlus
      chr.ccoords.ranLocPlus <- chr.ccoords.ranLocPlus+ccoords.ranLocPlus
      print(j)
    }
    a.coords.targetPlus <- a.coords.targetPlus+chr.acoords.targetPlus
    t.coords.targetPlus <- t.coords.targetPlus+chr.tcoords.targetPlus
    g.coords.targetPlus <- g.coords.targetPlus+chr.gcoords.targetPlus
    c.coords.targetPlus <- c.coords.targetPlus+chr.ccoords.targetPlus
    print(i)
    a.coords.ranLocPlus <- a.coords.ranLocPlus+chr.acoords.ranLocPlus
    t.coords.ranLocPlus <- t.coords.ranLocPlus+chr.tcoords.ranLocPlus
    g.coords.ranLocPlus <- g.coords.ranLocPlus+chr.gcoords.ranLocPlus
    c.coords.ranLocPlus <- c.coords.ranLocPlus+chr.ccoords.ranLocPlus
    print(i)

    # "Minus" strand loci
    chr.acoords.targetMinus <- rep(0, times = locSize)
    chr.tcoords.targetMinus <- rep(0, times = locSize)
    chr.gcoords.targetMinus <- rep(0, times = locSize)
    chr.ccoords.targetMinus <- rep(0, times = locSize)
    chr.acoords.ranLocMinus <- rep(0, times = locSize)
    chr.tcoords.ranLocMinus <- rep(0, times = locSize)
    chr.gcoords.ranLocMinus <- rep(0, times = locSize)
    chr.ccoords.ranLocMinus <- rep(0, times = locSize)
    print(i)
    for(j in 1:length(targetMinusStart)) {
      acoords.targetMinus <- rep(0, times = locSize)
      tcoords.targetMinus <- rep(0, times = locSize)
      gcoords.targetMinus <- rep(0, times = locSize)
      ccoords.targetMinus <- rep(0, times = locSize)
      print(j)

      sel.seq.targetMinus <- unlist(strsplit(as.character(chr.seq[targetMinusStart[j]:targetMinusEnd[j]]), split=""))
      acoords.targetMinus[which(sel.seq.targetMinus == "A")] <- 1
      tcoords.targetMinus[which(sel.seq.targetMinus == "T")] <- 1
      gcoords.targetMinus[which(sel.seq.targetMinus == "G")] <- 1
      ccoords.targetMinus[which(sel.seq.targetMinus == "C")] <- 1
      chr.acoords.targetMinus <- chr.acoords.targetMinus+acoords.targetMinus
      chr.tcoords.targetMinus <- chr.tcoords.targetMinus+tcoords.targetMinus
      chr.gcoords.targetMinus <- chr.gcoords.targetMinus+gcoords.targetMinus
      chr.ccoords.targetMinus <- chr.ccoords.targetMinus+ccoords.targetMinus
      print(j)

      acoords.ranLocMinus <- rep(0, times = locSize)
      tcoords.ranLocMinus <- rep(0, times = locSize)
      gcoords.ranLocMinus <- rep(0, times = locSize)
      ccoords.ranLocMinus <- rep(0, times = locSize)
      print(j)

      sel.seq.ranLocMinus <- unlist(strsplit(as.character(chr.seq[ranLocMinusStart[j]:ranLocMinusEnd[j]]), split=""))
      acoords.ranLocMinus[which(sel.seq.ranLocMinus == "A")] <- 1
      tcoords.ranLocMinus[which(sel.seq.ranLocMinus == "T")] <- 1
      gcoords.ranLocMinus[which(sel.seq.ranLocMinus == "G")] <- 1
      ccoords.ranLocMinus[which(sel.seq.ranLocMinus == "C")] <- 1
      chr.acoords.ranLocMinus <- chr.acoords.ranLocMinus+acoords.ranLocMinus
      chr.tcoords.ranLocMinus <- chr.tcoords.ranLocMinus+tcoords.ranLocMinus
      chr.gcoords.ranLocMinus <- chr.gcoords.ranLocMinus+gcoords.ranLocMinus
      chr.ccoords.ranLocMinus <- chr.ccoords.ranLocMinus+ccoords.ranLocMinus
      print(j)
    }
    a.coords.targetMinus <- a.coords.targetMinus+chr.acoords.targetMinus
    t.coords.targetMinus <- t.coords.targetMinus+chr.tcoords.targetMinus
    g.coords.targetMinus <- g.coords.targetMinus+chr.gcoords.targetMinus
    c.coords.targetMinus <- c.coords.targetMinus+chr.ccoords.targetMinus
    print(i)
    a.coords.ranLocMinus <- a.coords.ranLocMinus+chr.acoords.ranLocMinus
    t.coords.ranLocMinus <- t.coords.ranLocMinus+chr.tcoords.ranLocMinus
    g.coords.ranLocMinus <- g.coords.ranLocMinus+chr.gcoords.ranLocMinus
    c.coords.ranLocMinus <- c.coords.ranLocMinus+chr.ccoords.ranLocMinus
    print(i)
  }
  
  # "Plus" strand loci
  targetPlus.tot <- sum(targetPlus.chr.tots)
  print(targetPlus.tot)
  ranLocPlus.tot <- sum(ranLocPlus.chr.tots)
  print(ranLocPlus.tot)

  a.coords.targetPlus <- a.coords.targetPlus/targetPlus.tot
  t.coords.targetPlus <- t.coords.targetPlus/targetPlus.tot
  g.coords.targetPlus <- g.coords.targetPlus/targetPlus.tot
  c.coords.targetPlus <- c.coords.targetPlus/targetPlus.tot
  a.coords.ranLocPlus <- a.coords.ranLocPlus/ranLocPlus.tot
  t.coords.ranLocPlus <- t.coords.ranLocPlus/ranLocPlus.tot
  g.coords.ranLocPlus <- g.coords.ranLocPlus/ranLocPlus.tot
  c.coords.ranLocPlus <- c.coords.ranLocPlus/ranLocPlus.tot

  write.table(a.coords.targetPlus, file = paste0(outDir, targetName, "_targetPlus.a.txt"))
  write.table(t.coords.targetPlus, file = paste0(outDir, targetName, "_targetPlus.t.txt"))
  write.table(g.coords.targetPlus, file = paste0(outDir, targetName, "_targetPlus.g.txt"))
  write.table(c.coords.targetPlus, file = paste0(outDir, targetName, "_targetPlus.c.txt"))
  write.table(a.coords.ranLocPlus, file = paste0(outDir, targetName, "_ranLocPlus.a.txt"))
  write.table(t.coords.ranLocPlus, file = paste0(outDir, targetName, "_ranLocPlus.t.txt"))
  write.table(g.coords.ranLocPlus, file = paste0(outDir, targetName, "_ranLocPlus.g.txt"))
  write.table(c.coords.ranLocPlus, file = paste0(outDir, targetName, "_ranLocPlus.c.txt"))

  # "Minus" strand loci
  targetMinus.tot <- sum(targetMinus.chr.tots)
  print(targetMinus.tot)
  ranLocMinus.tot <- sum(ranLocMinus.chr.tots)
  print(ranLocMinus.tot)

  a.coords.targetMinus <- a.coords.targetMinus/targetMinus.tot
  t.coords.targetMinus <- t.coords.targetMinus/targetMinus.tot
  g.coords.targetMinus <- g.coords.targetMinus/targetMinus.tot
  c.coords.targetMinus <- c.coords.targetMinus/targetMinus.tot
  a.coords.ranLocMinus <- a.coords.ranLocMinus/ranLocMinus.tot
  t.coords.ranLocMinus <- t.coords.ranLocMinus/ranLocMinus.tot
  g.coords.ranLocMinus <- g.coords.ranLocMinus/ranLocMinus.tot
  c.coords.ranLocMinus <- c.coords.ranLocMinus/ranLocMinus.tot

  write.table(rev(a.coords.targetMinus), file = paste0(outDir, targetName, "_targetMinus.a.txt"))
  write.table(rev(t.coords.targetMinus), file = paste0(outDir, targetName, "_targetMinus.t.txt"))
  write.table(rev(g.coords.targetMinus), file = paste0(outDir, targetName, "_targetMinus.g.txt"))
  write.table(rev(c.coords.targetMinus), file = paste0(outDir, targetName, "_targetMinus.c.txt"))
  write.table(rev(a.coords.ranLocMinus), file = paste0(outDir, targetName, "_ranLocMinus.a.txt"))
  write.table(rev(t.coords.ranLocMinus), file = paste0(outDir, targetName, "_ranLocMinus.t.txt"))
  write.table(rev(g.coords.ranLocMinus), file = paste0(outDir, targetName, "_ranLocMinus.g.txt"))
  write.table(rev(c.coords.ranLocMinus), file = paste0(outDir, targetName, "_ranLocMinus.c.txt"))

  # Both strands
  a.coords.target <- (a.coords.targetPlus+rev(a.coords.targetMinus))/2
  t.coords.target <- (t.coords.targetPlus+rev(t.coords.targetMinus))/2
  g.coords.target <- (g.coords.targetPlus+rev(g.coords.targetMinus))/2
  c.coords.target <- (c.coords.targetPlus+rev(c.coords.targetMinus))/2
  a.coords.ranLoc <- (a.coords.ranLocPlus+rev(a.coords.ranLocMinus))/2
  t.coords.ranLoc <- (t.coords.ranLocPlus+rev(t.coords.ranLocMinus))/2
  g.coords.ranLoc <- (g.coords.ranLocPlus+rev(g.coords.ranLocMinus))/2
  c.coords.ranLoc <- (c.coords.ranLocPlus+rev(c.coords.ranLocMinus))/2

  write.table(a.coords.target, file = paste0(outDir, targetName, "_target.a.txt"))
  write.table(t.coords.target, file = paste0(outDir, targetName, "_target.t.txt"))
  write.table(g.coords.target, file = paste0(outDir, targetName, "_target.g.txt"))
  write.table(c.coords.target, file = paste0(outDir, targetName, "_target.c.txt"))
  write.table(a.coords.ranLoc, file = paste0(outDir, targetName, "_ranLoc.a.txt"))
  write.table(t.coords.ranLoc, file = paste0(outDir, targetName, "_ranLoc.t.txt"))
  write.table(g.coords.ranLoc, file = paste0(outDir, targetName, "_ranLoc.g.txt"))
  write.table(c.coords.ranLoc, file = paste0(outDir, targetName, "_ranLoc.c.txt"))
}

