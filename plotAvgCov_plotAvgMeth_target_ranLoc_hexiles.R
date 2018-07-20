# Function to plot mean coverage profile of one chromatin mark around target and random loci grouped into hexiles
plotAvgCov_hexiles <- function(xplot, dat1, ranDat1, flankSize, winSize, Ylabel1, flankLabL, flankLabR, startLab1, endLab1, startLab2, endLab2, mycols) {
  # target loci
  plot(xplot, dat1[,1],
       ylim = c(min(dat1[,1], dat1[,2], dat1[,3], dat1[,4], dat1[,5], dat1[,6],
                    ranDat1[,1], ranDat1[,2], ranDat1[,3], ranDat1[,4], ranDat1[,5], ranDat1[,6]),
                max(dat1[,1], dat1[,2], dat1[,3], dat1[,4], dat1[,5], dat1[,6],
                    ranDat1[,1], ranDat1[,2], ranDat1[,3], ranDat1[,4], ranDat1[,5], ranDat1[,6])),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, dat1[,2], col=mycols[2], lwd = 1.5)
  lines(xplot, dat1[,3], col=mycols[3], lwd = 1.5)
  lines(xplot, dat1[,4], col=mycols[4], lwd = 1.5)
  lines(xplot, dat1[,5], col=mycols[5], lwd = 1.5)
  lines(xplot, dat1[,6], col=mycols[6], lwd = 1.5)
  axis(side = 2, at = pretty(c(dat1[,1], dat1[,2], dat1[,3], dat1[,4], dat1[,5], dat1[,6],
                               ranDat1[,1], ranDat1[,2], ranDat1[,3], ranDat1[,4], ranDat1[,5], ranDat1[,6])))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1)
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1[,1])-(flankSize/winSize), length(dat1[,1])), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1[,1])-(flankSize/winSize), length(dat1[,1])), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1[,1])-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)

  # random loci
  plot(xplot, ranDat1[,1],
       ylim = c(min(dat1[,1], dat1[,2], dat1[,3], dat1[,4], dat1[,5], dat1[,6],
                    ranDat1[,1], ranDat1[,2], ranDat1[,3], ranDat1[,4], ranDat1[,5], ranDat1[,6]),
                max(dat1[,1], dat1[,2], dat1[,3], dat1[,4], dat1[,5], dat1[,6],
                    ranDat1[,1], ranDat1[,2], ranDat1[,3], ranDat1[,4], ranDat1[,5], ranDat1[,6])),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, ranDat1[,2], col=mycols[2], lwd = 1.5)
  lines(xplot, ranDat1[,3], col=mycols[3], lwd = 1.5)
  lines(xplot, ranDat1[,4], col=mycols[4], lwd = 1.5)
  lines(xplot, ranDat1[,5], col=mycols[5], lwd = 1.5)
  lines(xplot, ranDat1[,6], col=mycols[6], lwd = 1.5)
  axis(side = 2, at = pretty(c(dat1[,1], dat1[,2], dat1[,3], dat1[,4], dat1[,5], dat1[,6],
                               ranDat1[,1], ranDat1[,2], ranDat1[,3], ranDat1[,4], ranDat1[,5], ranDat1[,6])))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1[,1])-(flankSize/winSize), length(ranDat1[,1])), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1[,1])-(flankSize/winSize), length(ranDat1[,1])), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1[,1])-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
}

