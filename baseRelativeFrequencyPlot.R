# Function to plot summed A+T AND G+C relative frequencies around target and random loci 
mergeBaseFreqPlot <- function(at.coords, gc.coords, at.ran.coords, gc.ran.coords, flankSize, flankLabL, flankLabR, midpointLab1, midpointLab2, mycols, xplot, mainTitle1, mainTitle2) {
  # targets
  plot(xplot, at.coords, col = mycols[2], lwd = 1.5, type = "l",
       ylim = c(min(at.coords, gc.coords, at.ran.coords, gc.ran.coords),
                max(at.coords, gc.coords, at.ran.coords, gc.ran.coords)),
       ann = F, xaxt = "n")
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle1)
  mtext(side = 2, line = 2, cex = 0.8, text = "Base relative frequency")
  lines(xplot, gc.coords, col = mycols[1], lwd = 1.5)
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab1, flankLabR))
  axis(side = 2, at = pretty(c(at.coords, gc.coords, at.ran.coords, gc.ran.coords)))
  abline(v = 0, lty = 3)
  legend("right",
         legend = c("A+T", "G+C"),
         col = c(mycols[2], mycols[1]),
         text.col = c(mycols[2], mycols[1]),
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
  # ranLoc
  plot(xplot, at.ran.coords, col = mycols[2], lwd = 1.5, type = "l",
       ylim = c(min(at.coords, gc.coords, at.ran.coords, gc.ran.coords),
                max(at.coords, gc.coords, at.ran.coords, gc.ran.coords)),
       ann = F, xaxt = "n")
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle2)
  lines(xplot, gc.ran.coords, col = mycols[1], lwd = 1.5)
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab2, flankLabR))
  axis(side = 2, at = pretty(c(at.coords, gc.coords, at.ran.coords, gc.ran.coords)))
  abline(v = 0, lty = 3)
  legend("right",
         legend = c("A+T", "G+C"),
         col = c(mycols[2], mycols[1]),
         text.col = c(mycols[2], mycols[1]),
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}

# Function to plot summed A+T AND G+C relative frequencies on DIFFERENT Y-AXES around target and random loci 
mergeBaseFreqPlotDiffY <- function(at.coords, gc.coords, at.ran.coords, gc.ran.coords, flankSize, flankLabL, flankLabR, midpointLab1, midpointLab2, mycols, xplot, mainTitle1, mainTitle2) {
  # targets
  plot(xplot, at.coords, col = mycols[2], lwd = 1.5, type = "l",
       ylim = c(min(at.coords, at.ran.coords),
                max(at.coords, at.ran.coords)),
       ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(at.coords, at.ran.coords)))
  mtext(side = 2, line = 2, cex = 0.8, text = "A+T relative frequency", col = mycols[2])
  par(new = T)
  plot(xplot, gc.coords, col = mycols[1], lwd = 1.5, type = "l",
       ylim = c(min(gc.coords, gc.ran.coords),
                max(gc.coords, gc.ran.coords)),
       ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(gc.coords, gc.ran.coords)))
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab1, flankLabR))
  abline(v = 0, lty = 3)
  box(lwd = 1.5)
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle1)

  # ranLoc
  plot(xplot, at.ran.coords, col = mycols[2], lwd = 1.5, type = "l",
       ylim = c(min(at.coords, at.ran.coords),
                max(at.coords, at.ran.coords)),
       ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(at.coords, at.ran.coords)))
  par(new = T)
  plot(xplot, gc.ran.coords, col = mycols[1], lwd = 1.5, type = "l",
       ylim = c(min(gc.coords, gc.ran.coords),
                max(gc.coords, gc.ran.coords)),
       ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(gc.coords, gc.ran.coords)))
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab2, flankLabR))
  mtext(side = 4, line = 2, cex = 0.8, text = "G+C relative frequency", col = mycols[1])
  abline(v = 0, lty = 3)
  box(lwd = 1.5)
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle2)
}

# Function to plot summed A+T OR G+C relative frequencies around target and random loci 
ATorGCmergeBaseFreqPlot <- function(coords, ran.coords, flankSize, flankLabL, flankLabR, midpointLab1, midpointLab2, mycols, xplot, mainTitle1, mainTitle2, legendLab) {
  # targets
  plot(xplot, coords, col = mycols, lwd = 1.5, type = "l",
       ylim = c(min(coords, ran.coords),
                max(coords, ran.coords)),
       ann = F, xaxt = "n")
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle1)
  mtext(side = 2, line = 2, cex = 0.8, text = "Base relative frequency")
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab1, flankLabR))
  axis(side = 2, at = pretty(c(coords, ran.coords)))
  abline(v = 0, lty = 3)
  legend("right",
         legend = legendLab,
         col = c(mycols),
         text.col = c(mycols),
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
  # ranLoc
  plot(xplot, ran.coords, col = mycols, lwd = 1.5, type = "l",
       ylim = c(min(coords, ran.coords),
                max(coords, ran.coords)),
       ann = F, xaxt = "n")
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle2)
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab2, flankLabR))
  axis(side = 2, at = pretty(c(coords, ran.coords)))
  abline(v = 0, lty = 3)
  legend("right",
         legend = legendLab,
         col = c(mycols),
         text.col = c(mycols),
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}

# Function to plot base relative frequencies around target and random loci
baseFreqPlot <- function(coords, ran.coords, flankSize, flankLabL, flankLabR, midpointLab1, midpointLab2, mycols, xplot, mainTitle1, mainTitle2) {
  # targets
  plot(xplot, coords[[1]], col = mycols[1], lwd = 1.5, type = "l",
       ylim = c(min(coords[[1]], coords[[2]], coords[[3]], coords[[4]], ran.coords[[1]], ran.coords[[2]], ran.coords[[3]], ran.coords[[4]]),
                max(coords[[1]], coords[[2]], coords[[3]], coords[[4]], ran.coords[[1]], ran.coords[[2]], ran.coords[[3]], ran.coords[[4]])),
       ann = F, xaxt = "n", yaxt = "n")
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle1)
  mtext(side = 2, line = 2, cex = 0.8, text = "Base relative frequency")
  lines(xplot, coords[[2]], col = mycols[2], lwd = 1.5)
  lines(xplot, coords[[3]], col = mycols[3], lwd = 1.5)
  lines(xplot, coords[[4]], col = mycols[4], lwd = 1.5)
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab1, flankLabR))
  axis(side = 2, at = pretty(c(coords[[1]], coords[[2]], coords[[3]], coords[[4]], ran.coords[[1]], ran.coords[[2]], ran.coords[[3]], ran.coords[[4]])))
  abline(v = 0, lty = 3)
  legend("right",
         legend = c("A", "T", "G", "C"),
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
  # ranLoc
  plot(xplot, ran.coords[[1]], col = mycols[1], lwd = 1.5, type = "l",
       ylim = c(min(coords[[1]], coords[[2]], coords[[3]], coords[[4]], ran.coords[[1]], ran.coords[[2]], ran.coords[[3]], ran.coords[[4]]),
                max(coords[[1]], coords[[2]], coords[[3]], coords[[4]], ran.coords[[1]], ran.coords[[2]], ran.coords[[3]], ran.coords[[4]])),
       ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, ran.coords[[2]], col = mycols[2], lwd = 1.5)
  lines(xplot, ran.coords[[3]], col = mycols[3], lwd = 1.5)
  lines(xplot, ran.coords[[4]], col = mycols[4], lwd = 1.5)
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle2)
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab2, flankLabR))
  axis(side = 2, at = pretty(c(coords[[1]], coords[[2]], coords[[3]], coords[[4]], ran.coords[[1]], ran.coords[[2]], ran.coords[[3]], ran.coords[[4]])))
  abline(v = 0, lty = 3)
  legend("right",
         legend = c("A", "T", "G", "C"),
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}

# Function to plot base relative frequencies on DIFFERENT Y-AXES (A or T AND G or C) around target and random loci
baseFreqPlotDiffY <- function(coords, ran.coords, flankSize, flankLabL, flankLabR, midpointLab1, midpointLab2, mycols, xplot, mainTitle1, mainTitle2) {
  # targets
  plot(xplot, coords[[1]], col = mycols[1], lwd = 1.5, type = "l",
       ylim = c(min(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]]),
                max(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]])),
       ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, coords[[2]], col = mycols[2], lwd = 1.5)
  axis(side = 2, at = pretty(c(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]])))
  mtext(side = 2, line = 2, cex = 0.8, text = "A or T relative frequency", col = "blue")
  par(new = T)
  plot(xplot, coords[[3]], col = mycols[3], lwd = 1.5, type = "l",
       ylim = c(min(coords[[3]], coords[[4]], ran.coords[[3]], ran.coords[[4]]),
                max(coords[[3]], coords[[4]], ran.coords[[3]], ran.coords[[4]])),
       ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, coords[[4]], col = mycols[4], lwd = 1.5)
  axis(side = 4, at = pretty(c(coords[[3]], coords[[4]], ran.coords[[3]], ran.coords[[4]])))
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab1, flankLabR))
  abline(v = 0, lty = 3)
  legend("right",
         legend = c("A", "T", "G", "C"),
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle1)

  # ranLoc
  plot(xplot, ran.coords[[1]], col = mycols[1], lwd = 1.5, type = "l",
       ylim = c(min(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]]),
                max(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]])),
       ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, ran.coords[[2]], col = mycols[2], lwd = 1.5)
  axis(side = 2, at = pretty(c(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]])))
  par(new = T)
  plot(xplot, ran.coords[[3]], col = mycols[3], lwd = 1.5, type = "l",
       ylim = c(min(coords[[3]], coords[[4]], ran.coords[[3]], ran.coords[[4]]),
                max(coords[[3]], coords[[4]], ran.coords[[3]], ran.coords[[4]])),
       ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, ran.coords[[4]], col = mycols[4], lwd = 1.5)
  axis(side = 4, at = pretty(c(coords[[3]], coords[[4]], ran.coords[[3]], ran.coords[[4]])))
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab2, flankLabR))
  mtext(side = 4, line = 2, cex = 0.8, text = "G or C relative frequency", col = "red")
  abline(v = 0, lty = 3)
  legend("right",
         legend = c("A", "T", "G", "C"),
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle2)
}

# Function to plot separate A and T OR G and C base relative frequencies around target and random loci
ATorGCbaseFreqPlot <- function(coords, ran.coords, flankSize, flankLabL, flankLabR, midpointLab1, midpointLab2, mycols, xplot, mainTitle1, mainTitle2, legendLab) {
  # targets
  plot(xplot, coords[[1]], col = mycols[1], lwd = 1.5, type = "l",
       ylim = c(min(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]]),
                max(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]])),
       ann = F, xaxt = "n")
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle1)
  mtext(side = 2, line = 2, cex = 0.8, text = "Base relative frequency")
  lines(xplot, coords[[2]], col = mycols[2], lwd = 1.5)
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab1, flankLabR))
  axis(side = 2, at = pretty(c(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]])))
  abline(v = 0, lty = 3)
  legend("right",
         legend = legendLab,
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
  # ranLoc
  plot(xplot, ran.coords[[1]], col = mycols[1], lwd = 1.5, type = "l",
       ylim = c(min(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]]),
                max(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]])),
       ann = F, xaxt = "n")
  lines(xplot, ran.coords[[2]], col = mycols[2], lwd = 1.5)
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle2)
  axis(side = 1, at = c(-flankSize, 0, flankSize), labels = c(flankLabL, midpointLab2, flankLabR))
  axis(side = 2, at = pretty(c(coords[[1]], coords[[2]], ran.coords[[1]], ran.coords[[2]])))
  abline(v = 0, lty = 3)
  legend("right",
         legend = legendLab,
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}


