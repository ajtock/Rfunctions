# Function to plot mean coverage profile of one chromatin mark vs another around
# target and random loci (single-nucleotide loci)
plotAvgCovSN_oneVanother <- function(xplot, dat1, dat2, ranDat1, ranDat2,
                                     locTitle, ranLocTitle,
                                     flankSize, winSize, Ylabel,
                                     flankLabL, flankLabR,
                                     startLab1, startLab2,
                                     legendLoc, legendLab,
                                     mycols) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, dat2, ranDat1, ranDat2),
                max(dat1, dat2, ranDat1, ranDat2)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  mtext(side = 3, line = 0.1, text = locTitle)
  axis(side = 2, at = pretty(c(dat1, dat2, ranDat1, ranDat2)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel, col = mycols[1])
  lines(xplot, dat2, type = "l", lwd = 1.5, col = mycols[2])
  axis(side = 1, at = c(1, (flankSize/winSize)+0.5, length(dat1)), labels = c("", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+0.5, length(dat1)),
        text = c(flankLabL, startLab1, flankLabR))
  abline(v = (flankSize/winSize)+0.5, lty = 3)
  legend(legendLoc,
         legend = legendLab,
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
  box(lwd = 1.5)

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, dat2, ranDat1, ranDat2),
                max(dat1, dat2, ranDat1, ranDat2)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  mtext(side = 3, line = 0.1, text = ranLocTitle)
  axis(side = 2, at = pretty(c(dat1, dat2, ranDat1, ranDat2)))
  lines(xplot, ranDat2, type = "l", lwd = 1.5, col = mycols[2])
  axis(side = 1, at = c(1, (flankSize/winSize)+0.5, length(ranDat1)), labels = c("", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+0.5, length(ranDat1)),
        text = c(flankLabL, startLab2, flankLabR))
  abline(v = (flankSize/winSize)+0.5, lty = 3)
  legend(legendLoc,
         legend = legendLab,
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}

