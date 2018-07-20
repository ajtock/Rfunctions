# Function to plot mean coverage profile of one chromatin mark vs another using same
# Y-axis around target and random loci
plotAvgCov_1v2_1Y <- function(xplot, dat1, dat2, ranDat1, ranDat2,
                              locTitle, ranLocTitle,
                              flankSize, winSize, Ylabel,
                              flankLabL, flankLabR,
                              startLab1, endLab1,
                              startLab2, endLab2,
                              legendLoc, legendLab,
                              mycols) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, dat2, ranDat1, ranDat2),
                max(dat1, dat2, ranDat1, ranDat2)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  mtext(side = 3, line = 0.1, text = locTitle)
  axis(side = 2, at = pretty(c(dat1, dat2, ranDat1, ranDat2)), cex = 0.7)
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel, col = mycols[1])
  lines(xplot, dat2, type = "l", lwd = 1.5, col = mycols[2])
  axis(side = 1, at = c(1,
                        (flankSize/winSize)+1,
                        length(dat1)-(flankSize/winSize),
                        length(dat1)),
       labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1,
                                              (flankSize/winSize)+1,
                                              length(dat1)-(flankSize/winSize),
                                              length(dat1)),
        text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
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
  axis(side = 2, at = pretty(c(dat1, dat2, ranDat1, ranDat2)), cex = 0.7)
  lines(xplot, ranDat2, type = "l", lwd = 1.5, col = mycols[2])
  axis(side = 1, at = c(1,
                        (flankSize/winSize)+1,
                        length(ranDat1)-(flankSize/winSize),
                        length(ranDat1)),
       labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1,
                                              (flankSize/winSize)+1,
                                              length(ranDat1)-(flankSize/winSize),
                                              length(ranDat1)),
        text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  legend(legendLoc,
         legend = legendLab,
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}

