# Function to plot mean coverage profile of one chromatin mark vs another around target and random loci
# Same Y axis using dat1 and dat3, and dat2 and dat4
plotAvgCov_oneVanother_sameY <- function(xplot,
                                         dat1, dat2,
                                         ranDat1, ranDat2,
                                         dat3, dat4,
                                         ranDat3, ranDat4,
                                         flankSize, winSize,
                                         Ylabel1, Ylabel2,
                                         flankLabL, flankLabR,
                                         startLab1, endLab1,
                                         startLab2, endLab2,
                                         mycols) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1,
                    dat3, ranDat3),
                max(dat1, ranDat1,
                    dat3, ranDat3)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1, ranDat1,
                               dat3, ranDat3)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycols[1])
  par(new = T)
  plot(xplot, dat2,
       ylim = c(min(dat2, ranDat2,
                    dat4, ranDat4),
                max(dat2, ranDat2,
                    dat4, ranDat4)),
       type = "l", lwd = 1.5, col = mycols[2], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2, ranDat2,
                               dat4, ranDat4)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, ranDat1,
                    dat3, ranDat3),
                max(dat1, ranDat1,
                    dat3, ranDat3)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1, ranDat1,
                               dat3, ranDat3)))
  par(new = T)
  plot(xplot, ranDat2,
       ylim = c(min(dat2, ranDat2,
                    dat4, ranDat4),
                max(dat2, ranDat2,
                    dat4, ranDat4)),
       type = "l", lwd = 1.5, col = mycols[2], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2, ranDat2,
                               dat4, ranDat4)))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycols[2])
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
}

# Function to plot mean coverage profile of one chromatin mark vs mean DNA methylation profile (all contexts) target and random loci
# Same Y axis using dat1 and dat3, and dat2 and dat4
plotAvgCov_plotAvgMeth_sameY <- function(xplot,
                                         dat1, CGmethDat1, CHGmethDat1, CHHmethDat1,
                                         ranDat1, CGmethRanDat1, CHGmethRanDat1, CHHmethRanDat1,
                                         dat2, CGmethDat2, CHGmethDat2, CHHmethDat2,
                                         ranDat2, CGmethRanDat2, CHGmethRanDat2, CHHmethRanDat2,
                                         flankSize, winSize,
                                         Ylabel1, Ylabel2,
                                         flankLabL, flankLabR,
                                         startLab1, endLab1,
                                         startLab2, endLab2,
                                         legendLoc,
                                         mycols, mycolsMeth) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1,
                    dat2, ranDat2),
                max(dat1, ranDat1,
                    dat2, ranDat2)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1, ranDat1,
                               dat2, ranDat2)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycols[1])
  par(new = T)
  plot(xplot, CGmethDat1,
       ylim = c(min(CGmethDat1, CHGmethDat1, CHHmethDat1,
                    CGmethRanDat1, CHGmethRanDat1, CHHmethRanDat1,
                    CGmethDat2, CHGmethDat2, CHHmethDat2,
                    CGmethRanDat2, CHGmethRanDat2, CHHmethRanDat2),
                max(CGmethDat1, CHGmethDat1, CHHmethDat1,
                    CGmethRanDat1, CHGmethRanDat1, CHHmethRanDat1,
                    CGmethDat2, CHGmethDat2, CHHmethDat2,
                    CGmethRanDat2, CHGmethRanDat2, CHHmethRanDat2)),
       type = "l", lwd = 1.5, col = mycolsMeth[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, CHGmethDat1, type = "l", lwd = 1.5, col = mycolsMeth[2])
  lines(xplot, CHHmethDat1, type = "l", lwd = 1.5, col = mycolsMeth[3])
  axis(side = 4, at = pretty(c(CGmethDat1, CHGmethDat1, CHHmethDat1,
                               CGmethRanDat1, CHGmethRanDat1, CHHmethRanDat1,
                               CGmethDat2, CHGmethDat2, CHHmethDat2,
                               CGmethRanDat2, CHGmethRanDat2, CHHmethRanDat2)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = c("CG", "CHG", "CHH"),
         col = mycolsMeth,
         text.col = mycolsMeth,
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, ranDat1,
                    dat2, ranDat2),
                max(dat1, ranDat1,
                    dat2, ranDat2)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1, ranDat1,
                               dat2, ranDat2)))
  par(new = T)
  plot(xplot, CGmethRanDat1,
       ylim = c(min(CGmethDat1, CHGmethDat1, CHHmethDat1,
                    CGmethRanDat1, CHGmethRanDat1, CHHmethRanDat1,
                    CGmethDat2, CHGmethDat2, CHHmethDat2,
                    CGmethRanDat2, CHGmethRanDat2, CHHmethRanDat2),
                max(CGmethDat1, CHGmethDat1, CHHmethDat1,
                    CGmethRanDat1, CHGmethRanDat1, CHHmethRanDat1,
                    CGmethDat2, CHGmethDat2, CHHmethDat2,
                    CGmethRanDat2, CHGmethRanDat2, CHHmethRanDat2)),
       type = "l", lwd = 1.5, col = mycolsMeth[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, CHGmethRanDat1, type = "l", lwd = 1.5, col = mycolsMeth[2])
  lines(xplot, CHHmethRanDat1, type = "l", lwd = 1.5, col = mycolsMeth[3])
  axis(side = 4, at = pretty(c(CGmethDat1, CHGmethDat1, CHHmethDat1,
                               CGmethRanDat1, CHGmethRanDat1, CHHmethRanDat1)))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycolsMeth[2])
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = c("CG", "CHG", "CHH"),
         col = mycolsMeth,
         text.col = mycolsMeth,
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

