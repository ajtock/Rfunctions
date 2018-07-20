# Function to plot mean coverage profile of one chromatin mark vs another around target and random loci
plotAvgCov_oneVanother <- function(xplot, dat1, dat2, ranDat1, ranDat2, flankSize, winSize, Ylabel1, Ylabel2, flankLabL, flankLabR, startLab1, endLab1, startLab2, endLab2, mycols) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1),
                max(dat1, ranDat1)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1, ranDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycols[1])
  par(new = T)
  plot(xplot, dat2,
       ylim = c(min(dat2, ranDat2),
                max(dat2, ranDat2)),
       type = "l", lwd = 1.5, col = mycols[2], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2, ranDat2)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, ranDat1),
                max(dat1, ranDat1)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1, ranDat1)))
  par(new = T)
  plot(xplot, ranDat2,
       ylim = c(min(dat2, ranDat2),
                max(dat2, ranDat2)),
       type = "l", lwd = 1.5, col = mycols[2], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2, ranDat2)))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycols[2])
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
}

# Function to plot mean coverage profile of one chromatin mark vs another (other in both wt and mutant on same Y) around target and random loci
plotAvgCov_oneWTVanotherWTmutant <- function(xplot, dat1, dat2, mutantDat2, ranDat1, ranDat2, mutantRanDat2, flankSize, winSize, Ylabel1, Ylabel2, flankLabL, flankLabR, startLab1, endLab1, startLab2, endLab2, legendLoc, legendLabs, mycolsDat1, mycolsDat2) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1),
                max(dat1, ranDat1)),
       type = "l", lwd = 1.5, col = mycolsDat1[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1, ranDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycolsDat1[1])
  par(new = T)
  plot(xplot, dat2,
       ylim = c(min(dat2, ranDat2, mutantDat2, mutantRanDat2),
                max(dat2, ranDat2, mutantDat2, mutantRanDat2)),
       type = "l", lwd = 1.5, col = mycolsDat2[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantDat2, type = "l", lwd = 1.5, col = mycolsDat2[2])
  axis(side = 4, at = pretty(c(dat2, ranDat2, mutantDat2, mutantRanDat2)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = mycolsDat2,
         text.col = mycolsDat2,
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, ranDat1),
                max(dat1, ranDat1)),
       type = "l", lwd = 1.5, col = mycolsDat1[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1, ranDat1)))
  par(new = T)
  plot(xplot, ranDat2,
       ylim = c(min(dat2, ranDat2, mutantDat2, mutantRanDat2),
                max(dat2, ranDat2, mutantDat2, mutantRanDat2)),
       type = "l", lwd = 1.5, col = mycolsDat2[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantRanDat2, type = "l", lwd = 1.5, col = mycolsDat2[2])
  axis(side = 4, at = pretty(c(dat2, ranDat2, mutantDat2, mutantRanDat2)))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycolsDat2[1])
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = mycolsDat2,
         text.col = mycolsDat2,
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot mean coverage profile of one chromatin mark vs another (in both wt and mutant on same Y-axes) around target and random loci
plotAvgCov_oneWTmutantVanotherWTmutant <- function(xplot,
                                                   dat1, dat2,
                                                   mutantDat1, mutantDat2,
                                                   ranDat1, ranDat2,
                                                   mutantRanDat1, mutantRanDat2,
                                                   flankSize, winSize,
                                                   Ylabel1, Ylabel2,
                                                   flankLabL, flankLabR,
                                                   startLab1, endLab1,
                                                   startLab2, endLab2,
                                                   legendLoc, legendLabs,
                                                   mycolsDat1, mycolsDat2) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1, mutantDat1, mutantRanDat1),
                max(dat1, ranDat1, mutantDat1, mutantRanDat1)),
       type = "l", lwd = 1.5, col = mycolsDat1[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantDat1, type = "l", lwd = 1.5, col = mycolsDat1[2])
  axis(side = 2, at = pretty(c(dat1, ranDat1, mutantDat1, mutantRanDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycolsDat1[1])
  par(new = T)
  plot(xplot, dat2,
       ylim = c(min(dat2, ranDat2, mutantDat2, mutantRanDat2),
                max(dat2, ranDat2, mutantDat2, mutantRanDat2)),
       type = "l", lwd = 1.5, col = mycolsDat2[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantDat2, type = "l", lwd = 1.5, col = mycolsDat2[2])
  axis(side = 4, at = pretty(c(dat2, ranDat2, mutantDat2, mutantRanDat2)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(mycolsDat1, mycolsDat2),
         text.col = c(mycolsDat1, mycolsDat2),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, ranDat1, mutantDat1, mutantRanDat1),
                max(dat1, ranDat1, mutantDat1, mutantRanDat1)),
       type = "l", lwd = 1.5, col = mycolsDat1[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantRanDat1, type = "l", lwd = 1.5, col = mycolsDat1[2])
  axis(side = 2, at = pretty(c(dat1, ranDat1, mutantDat1, mutantRanDat1)))
  par(new = T)
  plot(xplot, ranDat2,
       ylim = c(min(dat2, ranDat2, mutantDat2, mutantRanDat2),
                max(dat2, ranDat2, mutantDat2, mutantRanDat2)),
       type = "l", lwd = 1.5, col = mycolsDat2[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantRanDat2, type = "l", lwd = 1.5, col = mycolsDat2[2])
  axis(side = 4, at = pretty(c(dat2, ranDat2, mutantDat2, mutantRanDat2)))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycolsDat2[1])
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(mycolsDat1, mycolsDat2),
         text.col = c(mycolsDat1, mycolsDat2),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot mean coverage profile of one chromatin mark (in both wt and mutant on same Y-axes) vs another (in wt only) around target and random loci
plotAvgCov_oneWTmutantVanotherWT <- function(xplot,
                                             dat1, dat2,
                                             mutantDat1, 
                                             ranDat1, ranDat2,
                                             mutantRanDat1,
                                             flankSize, winSize,
                                             Ylabel1, Ylabel2,
                                             flankLabL, flankLabR,
                                             startLab1, endLab1,
                                             startLab2, endLab2,
                                             legendLoc, legendLabs,
                                             mycolsDat1, mycolsDat2) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1, mutantDat1, mutantRanDat1),
                max(dat1, ranDat1, mutantDat1, mutantRanDat1)),
       type = "l", lwd = 1.5, col = mycolsDat1[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantDat1, type = "l", lwd = 1.5, col = mycolsDat1[2])
  axis(side = 2, at = pretty(c(dat1, ranDat1, mutantDat1, mutantRanDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycolsDat1[1])
  par(new = T)
  plot(xplot, dat2,
       ylim = c(min(dat2, ranDat2),
                max(dat2, ranDat2)),
       type = "l", lwd = 1.5, col = mycolsDat2[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2, ranDat2)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(mycolsDat1, mycolsDat2),
         text.col = c(mycolsDat1, mycolsDat2),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, ranDat1, mutantDat1, mutantRanDat1),
                max(dat1, ranDat1, mutantDat1, mutantRanDat1)),
       type = "l", lwd = 1.5, col = mycolsDat1[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantRanDat1, type = "l", lwd = 1.5, col = mycolsDat1[2])
  axis(side = 2, at = pretty(c(dat1, ranDat1, mutantDat1, mutantRanDat1)))
  par(new = T)
  plot(xplot, ranDat2,
       ylim = c(min(dat2, ranDat2),
                max(dat2, ranDat2)),
       type = "l", lwd = 1.5, col = mycolsDat2[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2, ranDat2)))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycolsDat2[1])
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(mycolsDat1, mycolsDat2),
         text.col = c(mycolsDat1, mycolsDat2),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}


# Function to plot mean coverage profile of one chromatin mark vs another (in both wt and mutant on same Y-axes) around target and random loci,
# on same Y-axis as equivalent data for other targets and random loci
plotAvgCov_oneWTmutantVanotherWTmutant_sameY  <- function(xplot,
                                                          dat1, dat2, dat3, dat4,
                                                          mutantDat1, mutantDat2, mutantDat3, mutantDat4,
                                                          ranDat1, ranDat2, ranDat3, ranDat4,
                                                          mutantRanDat1, mutantRanDat2, mutantRanDat3, mutantRanDat4,
                                                          flankSize, winSize,
                                                          Ylabel1, Ylabel2,
                                                          flankLabL, flankLabR,
                                                          startLab1, endLab1,
                                                          startLab2, endLab2,
                                                          legendLoc, legendLabs,
                                                          mycolsDat1, mycolsDat2) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3),
                max(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3)),
       type = "l", lwd = 1.5, col = mycolsDat1[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantDat1, type = "l", lwd = 1.5, col = mycolsDat1[2])
  axis(side = 2, at = pretty(c(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycolsDat1[1])
  par(new = T)
  plot(xplot, dat2,
       ylim = c(min(dat2, ranDat2, mutantDat2, mutantRanDat2, dat4, ranDat4, mutantDat4, mutantRanDat4),
                max(dat2, ranDat2, mutantDat2, mutantRanDat2, dat4, ranDat4, mutantDat4, mutantRanDat4)),
       type = "l", lwd = 1.5, col = mycolsDat2[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantDat2, type = "l", lwd = 1.5, col = mycolsDat2[2])
  axis(side = 4, at = pretty(c(dat2, ranDat2, mutantDat2, mutantRanDat2, dat4, ranDat4, mutantDat4, mutantRanDat4)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(mycolsDat1, mycolsDat2),
         text.col = c(mycolsDat1, mycolsDat2),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3),
                max(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3)),
       type = "l", lwd = 1.5, col = mycolsDat1[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantRanDat1, type = "l", lwd = 1.5, col = mycolsDat1[2])
  axis(side = 2, at = pretty(c(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3)))
  par(new = T)
  plot(xplot, ranDat2,
       ylim = c(min(dat2, ranDat2, mutantDat2, mutantRanDat2, dat4, ranDat4, mutantDat4, mutantRanDat4),
                max(dat2, ranDat2, mutantDat2, mutantRanDat2, dat4, ranDat4, mutantDat4, mutantRanDat4)),
       type = "l", lwd = 1.5, col = mycolsDat2[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantRanDat2, type = "l", lwd = 1.5, col = mycolsDat2[2])
  axis(side = 4, at = pretty(c(dat2, ranDat2, mutantDat2, mutantRanDat2, dat4, ranDat4, mutantDat4, mutantRanDat4)))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycolsDat2[1])
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(mycolsDat1, mycolsDat2),
         text.col = c(mycolsDat1, mycolsDat2),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot mean coverage profile of one chromatin mark (in both wt and mutant on same Y-axes) vs another (in wt only) around target and random loci
# on same Y-axis as equivalent data for other targets and random loci
plotAvgCov_oneWTmutantVanotherWT_sameY <- function(xplot,
                                                   dat1, dat2, dat3, dat4,
                                                   mutantDat1, mutantDat3,
                                                   ranDat1, ranDat2, ranDat3, ranDat4,
                                                   mutantRanDat1, mutantRanDat3,
                                                   flankSize, winSize,
                                                   Ylabel1, Ylabel2,
                                                   flankLabL, flankLabR,
                                                   startLab1, endLab1,
                                                   startLab2, endLab2,
                                                   legendLoc, legendLabs,
                                                   mycolsDat1, mycolsDat2) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3),
                max(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3)),
       type = "l", lwd = 1.5, col = mycolsDat1[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantDat1, type = "l", lwd = 1.5, col = mycolsDat1[2])
  axis(side = 2, at = pretty(c(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycolsDat1[1])
  par(new = T)
  plot(xplot, dat2,
       ylim = c(min(dat2, ranDat2, dat4, ranDat4),
                max(dat2, ranDat2, dat4, ranDat4)),
       type = "l", lwd = 1.5, col = mycolsDat2[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2, ranDat2, dat4, ranDat4)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(mycolsDat1, mycolsDat2),
         text.col = c(mycolsDat1, mycolsDat2),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3),
                max(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3)),
       type = "l", lwd = 1.5, col = mycolsDat1[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantRanDat1, type = "l", lwd = 1.5, col = mycolsDat1[2])
  axis(side = 2, at = pretty(c(dat1, ranDat1, mutantDat1, mutantRanDat1, dat3, ranDat3, mutantDat3, mutantRanDat3)))
  par(new = T)
  plot(xplot, ranDat2,
       ylim = c(min(dat2, ranDat2, dat4, ranDat4),
                max(dat2, ranDat2, dat4, ranDat4)),
       type = "l", lwd = 1.5, col = mycolsDat2[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2, ranDat2, dat4, ranDat4)))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycolsDat2[1])
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(mycolsDat1, mycolsDat2),
         text.col = c(mycolsDat1, mycolsDat2),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}


# Function to plot mean coverage profile of one chromatin mark vs mean DNA methylation profile (all contexts) target and random loci
plotAvgCov_plotAvgMeth <- function(xplot, dat1, CGmethDat, CHGmethDat, CHHmethDat, ranDat1, CGmethRanDat, CHGmethRanDat, CHHmethRanDat, flankSize, winSize, Ylabel1, Ylabel2, flankLabL, flankLabR, startLab1, endLab1, startLab2, endLab2, legendLoc, mycols, mycolsMeth) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1),
                max(dat1, ranDat1)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1, ranDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycols[1])
  par(new = T)
  plot(xplot, CGmethDat,
       ylim = c(min(CGmethDat, CHGmethDat, CHHmethDat,
                    CGmethRanDat, CHGmethRanDat, CHHmethRanDat),
                max(CGmethDat, CHGmethDat, CHHmethDat,
                    CGmethRanDat, CHGmethRanDat, CHHmethRanDat)),
       type = "l", lwd = 1.5, col = mycolsMeth[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, CHGmethDat, type = "l", lwd = 1.5, col = mycolsMeth[2])
  lines(xplot, CHHmethDat, type = "l", lwd = 1.5, col = mycolsMeth[3])
  axis(side = 4, at = pretty(c(CGmethDat, CHGmethDat, CHHmethDat,
                               CGmethRanDat, CHGmethRanDat, CHHmethRanDat)))
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
       ylim = c(min(dat1, ranDat1),
                max(dat1, ranDat1)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1, ranDat1)))
  par(new = T)
  plot(xplot, CGmethRanDat,
       ylim = c(min(CGmethDat, CHGmethDat, CHHmethDat,
                    CGmethRanDat, CHGmethRanDat, CHHmethRanDat),
                max(CGmethDat, CHGmethDat, CHHmethDat,
                    CGmethRanDat, CHGmethRanDat, CHHmethRanDat)),
       type = "l", lwd = 1.5, col = mycolsMeth[1], ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, CHGmethRanDat, type = "l", lwd = 1.5, col = mycolsMeth[2])
  lines(xplot, CHHmethRanDat, type = "l", lwd = 1.5, col = mycolsMeth[3])
  axis(side = 4, at = pretty(c(CGmethDat, CHGmethDat, CHHmethDat,
                               CGmethRanDat, CHGmethRanDat, CHHmethRanDat)))
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
