# Function to plot mean coverage profile of one chromatin mark vs another around target and random loci
plotAvgCov_oneVanother <- function(xplot, dat1, dat2, ranDat1, ranDat2, flankSize, winSize, Ylabel1, Ylabel2, flankLabL, flankLabR, startLab1, endLab1, startLab2, endLab2, mycols) {
  # target loci
  plot(xplot, dat1[,1],
       ylim = c(min(dat1[,1], ranDat1[,1]),
                max(dat1[,1], ranDat1[,1])),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1[,1], ranDat1[,1])))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycols[1])
  par(new = T)
  plot(xplot, dat2[,1],
       ylim = c(min(dat2[,1], ranDat2[,1]),
                max(dat2[,1], ranDat2[,1])),
       type = "l", lwd = 1.5, col = mycols[2], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2[,1], ranDat2[,1])))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1[,1])-(flankSize/winSize), length(dat1[,1])), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1[,1])-(flankSize/winSize), length(dat1[,1])), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1[,1])-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)

  # random loci
  plot(xplot, ranDat1[,1],
       ylim = c(min(dat1[,1], ranDat1[,1]),
                max(dat1[,1], ranDat1[,1])),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1[,1], ranDat1[,1])))
  par(new = T)
  plot(xplot, ranDat2[,1],
       ylim = c(min(dat2[,1], ranDat2[,1]),
                max(dat2[,1], ranDat2[,1])),
       type = "l", lwd = 1.5, col = mycols[2], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2[,1], ranDat2[,1])))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycols[2])
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1[,1])-(flankSize/winSize), length(ranDat1[,1])), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1[,1])-(flankSize/winSize), length(ranDat1[,1])), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1[,1])-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
}

# Function to plot mean coverage profile of one chromatin mark vs mean DNA methylation profile (all contexts) target and random loci
plotAvgCov_plotAvgMeth <- function(xplot, dat1, CGmethDat, CHGmethDat, CHHmethDat, ranDat1, CGmethRanDat, CHGmethRanDat, CHHmethRanDat, flankSize, winSize, Ylabel1, Ylabel2, flankLabL, flankLabR, startLab1, endLab1, startLab2, endLab2, legendLoc, mycols, mycolsMeth) {
  # target loci
  plot(xplot, dat1[,1],
       ylim = c(min(dat1[,1], ranDat1[,1]),
                max(dat1[,1], ranDat1[,1])),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1[,1], ranDat1[,1])))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycols[1])
  par(new = T)
  plot(xplot, CGmethDat[,1],
       ylim = c(min(CGmethDat[,1], CHGmethDat[,1], CHHmethDat[,1],
                    CGmethRanDat[,1], CHGmethRanDat[,1], CHHmethRanDat[,1]),
                max(CGmethDat[,1], CHGmethDat[,1], CHHmethDat[,1],
                    CGmethRanDat[,1], CHGmethRanDat[,1], CHHmethRanDat[,1])),
       type = "l", lwd = 1.5, col = mycolsMeth[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, CHGmethDat[,1], type = "l", lwd = 1.5, col = mycolsMeth[2])
  lines(xplot, CHHmethDat[,1], type = "l", lwd = 1.5, col = mycolsMeth[3])
  axis(side = 4, at = pretty(c(CGmethDat[,1], CHGmethDat[,1], CHHmethDat[,1],
                               CGmethRanDat[,1], CHGmethRanDat[,1], CHHmethRanDat[,1])))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1[,1])-(flankSize/winSize), length(dat1[,1])), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1[,1])-(flankSize/winSize), length(dat1[,1])), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1[,1])-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = c("CG", "CHG", "CHH"),
         col = mycolsMeth,
         text.col = mycolsMeth,
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1[,1],
       ylim = c(min(dat1[,1], ranDat1[,1]),
                max(dat1[,1], ranDat1[,1])),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1[,1], ranDat1[,1])))
  par(new = T)
  plot(xplot, CGmethRanDat[,1],
       ylim = c(min(CGmethDat[,1], CHGmethDat[,1], CHHmethDat[,1],
                    CGmethRanDat[,1], CHGmethRanDat[,1], CHHmethRanDat[,1]),
                max(CGmethDat[,1], CHGmethDat[,1], CHHmethDat[,1],
                    CGmethRanDat[,1], CHGmethRanDat[,1], CHHmethRanDat[,1])),
       type = "l", lwd = 1.5, col = mycolsMeth[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, CHGmethRanDat[,1], type = "l", lwd = 1.5, col = mycolsMeth[2])
  lines(xplot, CHHmethRanDat[,1], type = "l", lwd = 1.5, col = mycolsMeth[3])
  axis(side = 4, at = pretty(c(CGmethDat[,1], CHGmethDat[,1], CHHmethDat[,1],
                               CGmethRanDat[,1], CHGmethRanDat[,1], CHHmethRanDat[,1])))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycolsMeth[2])
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1[,1])-(flankSize/winSize), length(ranDat1[,1])), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1[,1])-(flankSize/winSize), length(ranDat1[,1])), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1[,1])-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = c("CG", "CHG", "CHH"),
         col = mycolsMeth,
         text.col = mycolsMeth,
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}
