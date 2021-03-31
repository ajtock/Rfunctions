# Function to plot mean coverage profile of one chromatin mark vs another around target and random loci
plotAvgCov_oneVanother <- function(xplot, dat1, dat2, ranDat1, ranDat2, flankSize, winSize, Ylabel1, Ylabel2, flankLabL, flankLabR, startLab1, endLab1, startLab2, endLab2, mycols) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1),
                max(dat1, ranDat1)),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(dat1, ranDat1)))
  mtext(side = 2, line = 1.5, cex = 0.8, text = Ylabel1, col = mycols[1])
  par(new = T)
  plot(xplot, dat2,
       ylim = c(min(dat2, ranDat2),
                max(dat2, ranDat2)),
       type = "l", lwd = 1.5, col = mycols[2], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2, ranDat2)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 0.5, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
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
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -3.0), labels = Ylabel2, xpd = NA, srt = -90, col = mycols[2])
#  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycols[2])
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 0.5, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
}

# Function to plot mean coverage profile of one chromatin mark vs another around midpoint of targets and random loci (one Y-axis)
plotAvgCovMidpoint_1v2_oneY <- function(xplot,
                                        dat1, dat2,
                                        ranDat1, ranDat2,
                                        col1, col2,
                                        Ylabel,
                                        winSize,
                                        flankSize,
                                        flankLabL, flankLabR,
                                        legendLabs, legendLoc) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1,
                    dat2, ranDat2),
                max(dat1, ranDat1,
                    dat2, ranDat2)),
       type = "l", lwd = 3, col = col1,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s])*" = "*.(round(cor(dat1,
                                                      dat2,
                                                      method = "spearman"),
                                                      digits = 2))),
       cex.main = 1.5)
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel, col = "black")
  lines(xplot, dat2, type = "l", lwd = 3, col = col2)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5,
       at = c(1, (flankSize/winSize)+1, length(dat1)),
       labels = c("", "", ""))
  mtext(side = 1, line = 0.75, cex = 0.9,
        at = c(1, (flankSize/winSize)+1, length(dat1)),
        text = c(flankLabL, "Midpoint", flankLabR))
  abline(v = (flankSize/winSize)+1, lty = 3, lwd = 3)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2),
         text.col = c(col1, col2),
         text.font = c(1, 1),
         ncol = 1, cex = 0.4, lwd = 1.5, bty = "n")
  box(lwd = 1.5)

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, ranDat1,
                    dat2, ranDat2),
                max(dat1, ranDat1,
                    dat2, ranDat2)),
       type = "l", lwd = 3, col = col1,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s])*" = "*.(round(cor(ranDat1,
                                                      ranDat2,
                                                      method = "spearman"),
                                                      digits = 2))),
       cex.main = 1.5)
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel, col = "black")
  lines(xplot, ranDat2, type = "l", lwd = 3, col = col2)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5,
       at = c(1, (flankSize/winSize)+1, length(dat1)),
       labels = c("", "", ""))
  mtext(side = 1, line = 0.75, cex = 0.9,
        at = c(1, (flankSize/winSize)+1, length(dat1)),
        text = c(flankLabL, "Midpoint", flankLabR))
  abline(v = (flankSize/winSize)+1, lty = 3, lwd = 3)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2),
         text.col = c(col1, col2),
         text.font = c(1, 1),
         ncol = 1, cex = 0.4, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}

# Function to plot mean coverage profile of one chromatin mark vs another around midpoint of targets and random loci (one Y-axis)
plotAvgCovMidpoint_1v2_oneY_Ylim <- function(featureTitle,
                                             xplot,
                                             dat1, dat2,
                                             ranDat1, ranDat2,
                                             col1, col2,
                                             Ylim,
                                             Ylabel,
                                             winSize,
                                             flankSize,
                                             flankLabL, flankLabR,
                                             legendLabs,
                                             targetLegendLoc, ranLocLegendLoc) {
  # target loci
  plot(xplot, dat1,
       ylim = Ylim,
       type = "l", lwd = 3, col = col1,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(featureTitle)~~italic("r"[s])*" = "*.(round(cor(dat1,
                                                      dat2,
                                                      method = "spearman"),
                                                      digits = 2))),
       cex.main = 1.5)
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel, col = "black")
  lines(xplot, dat2, type = "l", lwd = 3, col = col2)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5,
       at = c(1, (flankSize/winSize)+1, length(dat1)),
       labels = c("", "", ""))
  mtext(side = 1, line = 0.75, cex = 0.9,
        at = c(1, (flankSize/winSize)+1, length(dat1)),
        text = c(flankLabL, "Midpoint", flankLabR))
  abline(v = (flankSize/winSize)+1, lty = 3, lwd = 3)
  legend(targetLegendLoc,
         legend = legendLabs,
         col = c(col1, col2),
         text.col = c(col1, col2),
         text.font = c(1, 1),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
  box(lwd = 1.5)

  # random loci
  plot(xplot, ranDat1,
       ylim = Ylim,
       type = "l", lwd = 3, col = col1,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote("Random loci"~~italic("r"[s])*" = "*.(round(cor(ranDat1,
                                                      ranDat2,
                                                      method = "spearman"),
                                                      digits = 2))),
       cex.main = 1.5)
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel, col = "black")
  lines(xplot, ranDat2, type = "l", lwd = 3, col = col2)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5,
       at = c(1, (flankSize/winSize)+1, length(dat1)),
       labels = c("", "", ""))
  mtext(side = 1, line = 0.75, cex = 0.9,
        at = c(1, (flankSize/winSize)+1, length(dat1)),
        text = c(flankLabL, "Midpoint", flankLabR))
  abline(v = (flankSize/winSize)+1, lty = 3, lwd = 3)
  legend(ranLocLegendLoc,
         legend = legendLabs,
         col = c(col1, col2),
         text.col = c(col1, col2),
         text.font = c(1, 1),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}

# Function to plot mean coverage profile of one chromatin mark vs another around midpoint of targets and random loci (two Y-axes)
plotAvgCovMidpoint_1v2_twoY_Ylims <- function(featureTitle,
                                              xplot,
                                              dat1, dat2,
                                              ranDat1, ranDat2,
                                              col1, col2,
                                              Ylim1, Ylim2,
                                              Ylabel1, Ylabel2,
                                              winSize,
                                              flankSize,
                                              flankLabL, flankLabR) {
  # target loci
  plot(xplot, dat1,
       ylim = Ylim1,
       type = "l", lwd = 3, col = col1,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(featureTitle)~~italic("r"[s])*" = "*.(round(cor(dat1,
                                                                       dat2,
                                                                       method = "spearman"),
                                                                       digits = 2))),
       cex.main = 1.5)
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = col1)
  par(new = T)
  plot(xplot, dat2,
       ylim = Ylim2,
       type = "l", lwd = 3, col = col2,
       ann = F,
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 1, lwd.tick = 1.5)
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       at = c(1, (flankSize/winSize)+1, length(dat1)),
       labels = c("", "", ""))
  mtext(side = 1, line = 0.75, cex = 0.9,
        at = c(1, (flankSize/winSize)+1, length(dat1)),
        text = c(flankLabL, "Midpoint", flankLabR))
  abline(v = (flankSize/winSize)+1, lty = 3, lwd = 3)
  box(lwd = 1.5)

  # random loci
  plot(xplot, ranDat1,
       ylim = Ylim1,
       type = "l", lwd = 3, col = col1,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote("Random loci"~~italic("r"[s])*" = "*.(round(cor(ranDat1,
                                                      ranDat2,
                                                      method = "spearman"),
                                                      digits = 2))),
       cex.main = 1.5)
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  par(new = T)
  plot(xplot, ranDat2,
       ylim = Ylim2,
       type = "l", lwd = 3, col = col2,
       ann = F,
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = col2)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5,
       at = c(1, (flankSize/winSize)+1, length(dat1)),
       labels = c("", "", ""))
  mtext(side = 1, line = 0.75, cex = 0.9,
        at = c(1, (flankSize/winSize)+1, length(dat1)),
        text = c(flankLabL, "Midpoint", flankLabR))
  abline(v = (flankSize/winSize)+1, lty = 3, lwd = 3)
  box(lwd = 1.5)
}

# Function to plot mean coverage profile of one chromatin mark vs another around midpoint of targets (one Y-axis)
plotAvgCovMidpoint_1v2_oneY_Ylim_noRan <- function(featureTitle,
                                                   xplot,
                                                   dat1, dat2,
                                                   col1, col2,
                                                   Ylim,
                                                   Ylabel,
                                                   winSize,
                                                   flankSize,
                                                   flankLabL, flankLabR,
                                                   legendLabs, legendLoc) {
  # target loci
  plot(xplot, dat1,
       ylim = Ylim,
       type = "l", lwd = 3, col = col1,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(featureTitle)~~italic("r"[s])*" = "*.(round(cor(dat1,
                                                      dat2,
                                                      method = "spearman"),
                                                      digits = 2))),
       cex.main = 1.5)
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel, col = "black")
  lines(xplot, dat2, type = "l", lwd = 3, col = col2)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5,
       at = c(1, (flankSize/winSize)+1, length(dat1)),
       labels = c("", "", ""))
  mtext(side = 1, line = 0.75, cex = 0.9,
        at = c(1, (flankSize/winSize)+1, length(dat1)),
        text = c(flankLabL, "Midpoint", flankLabR))
  abline(v = (flankSize/winSize)+1, lty = 3, lwd = 3)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2),
         text.col = c(col1, col2),
         text.font = c(1, 3),
         ncol = 1, cex = 1, lwd = 1.5, bty = "n")
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

# Function to plot mean coverage profile of wt vs mutant chromatin mark on same Y-axes around target and random loci
plotAvgCov_WTvMutant <- function(xplot,
                                 dat1,
                                 mutantDat1,
                                 ranDat1,
                                 mutantRanDat1,
                                 flankSize, winSize,
                                 Ylabel1,
                                 flankLabL, flankLabR,
                                 startLab1, endLab1,
                                 startLab2, endLab2,
                                 legendLoc, legendLabs,
                                 wtCol, mutantCol) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, ranDat1, mutantDat1, mutantRanDat1),
                max(dat1, ranDat1, mutantDat1, mutantRanDat1)),
       type = "l", lwd = 1.5, col = wtCol, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantDat1, type = "l", lwd = 1.5, col = mutantCol)
  axis(side = 2, at = pretty(c(dat1, ranDat1, mutantDat1, mutantRanDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = wtCol)
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(wtCol, mutantCol),
         text.col = c(wtCol, mutantCol),
         text.font = c(1, 3),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, ranDat1, mutantDat1, mutantRanDat1),
                max(dat1, ranDat1, mutantDat1, mutantRanDat1)),
       type = "l", lwd = 1.5, col = wtCol, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, mutantRanDat1, type = "l", lwd = 1.5, col = mutantCol)
  axis(side = 2, at = pretty(c(dat1, ranDat1, mutantDat1, mutantRanDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = wtCol)
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(wtCol, mutantCol),
         text.col = c(wtCol, mutantCol),
         text.font = c(1, 3),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot mean coverage profile of wt vs mutant chromatin mark for DNA and RNA TEs on same Y-axes around target and random loci
plotAvgCov_WTvMutant_DNA_RNA_TEs <- function(xplot,
                                             DNAdat1,
                                             DNAmutantDat1,
                                             DNAranDat1,
                                             DNAmutantRanDat1,
                                             RNAdat1,
                                             RNAmutantDat1,
                                             RNAranDat1,
                                             RNAmutantRanDat1,
                                             flankSize, winSize,
                                             Ylabel1,
                                             flankLabL, flankLabR,
                                             startLab1, endLab1,
                                             startLab2, endLab2,
                                             legendLoc, legendLabs,
                                             wtCol, mutantCol) {
  # target loci
  plot(xplot, DNAdat1,
       ylim = c(min(DNAdat1, DNAranDat1, DNAmutantDat1, DNAmutantRanDat1, RNAdat1, RNAranDat1, RNAmutantDat1, RNAmutantRanDat1),
                max(DNAdat1, DNAranDat1, DNAmutantDat1, DNAmutantRanDat1, RNAdat1, RNAranDat1, RNAmutantDat1, RNAmutantRanDat1)),
       type = "l", lwd = 1.5, col = wtCol, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, DNAmutantDat1, type = "l", lwd = 1.5, col = mutantCol)
  axis(side = 2, at = pretty(c(DNAdat1, DNAranDat1, DNAmutantDat1, DNAmutantRanDat1, RNAdat1, RNAranDat1, RNAmutantDat1, RNAmutantRanDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = wtCol)
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(DNAdat1)-(flankSize/winSize), length(DNAdat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(DNAdat1)-(flankSize/winSize), length(DNAdat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(DNAdat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(wtCol, mutantCol),
         text.col = c(wtCol, mutantCol),
         text.font = c(1, 3),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, DNAranDat1,
       ylim = c(min(DNAdat1, DNAranDat1, DNAmutantDat1, DNAmutantRanDat1, RNAdat1, RNAranDat1, RNAmutantDat1, RNAmutantRanDat1),
                max(DNAdat1, DNAranDat1, DNAmutantDat1, DNAmutantRanDat1, RNAdat1, RNAranDat1, RNAmutantDat1, RNAmutantRanDat1)),
       type = "l", lwd = 1.5, col = wtCol, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, DNAmutantRanDat1, type = "l", lwd = 1.5, col = mutantCol)
  axis(side = 2, at = pretty(c(DNAdat1, DNAranDat1, DNAmutantDat1, DNAmutantRanDat1, RNAdat1, RNAranDat1, RNAmutantDat1, RNAmutantRanDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = wtCol)
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(DNAranDat1)-(flankSize/winSize), length(DNAranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(DNAranDat1)-(flankSize/winSize), length(DNAranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(DNAranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(wtCol, mutantCol),
         text.col = c(wtCol, mutantCol),
         text.font = c(1, 3),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot mean coverage profile of one vs another chromatin mark for DNA and RNA TEs on separate Y-axes around target and random loci
plotAvgCov_oneVanother_DNA_RNA_TEs <- function(xplot,
                                               DNAdat1,
                                               DNAranDat1,
                                               RNAdat1,
                                               RNAranDat1,
                                               DNAdat2,
                                               DNAranDat2,
                                               RNAdat2,
                                               RNAranDat2,
                                               flankSize, winSize,
                                               Ylabel1, Ylabel2,
                                               flankLabL, flankLabR,
                                               startLab1, endLab1,
                                               startLab2, endLab2,
                                               dat1Col, dat2Col) {
  # target loci
  plot(xplot, DNAdat1,
       ylim = c(min(DNAdat1, DNAranDat1, RNAdat1, RNAranDat1),
                max(DNAdat1, DNAranDat1, RNAdat1, RNAranDat1)),
       type = "l", lwd = 1.5, col = dat1Col, ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(DNAdat1, DNAranDat1, RNAdat1, RNAranDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = dat1Col)
  par(new = T)
  plot(xplot, DNAdat2,
       ylim = c(min(DNAdat2, DNAranDat2, RNAdat2, RNAranDat2),
                max(DNAdat2, DNAranDat2, RNAdat2, RNAranDat2)),
       type = "l", lwd = 1.5, col = dat2Col, ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(DNAdat2, DNAranDat2, RNAdat2, RNAranDat2)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(DNAdat1)-(flankSize/winSize), length(DNAdat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(DNAdat1)-(flankSize/winSize), length(DNAdat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(DNAdat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)

  # random loci
  plot(xplot, DNAranDat1,
       ylim = c(min(DNAdat1, DNAranDat1, RNAdat1, RNAranDat1),
                max(DNAdat1, DNAranDat1, RNAdat1, RNAranDat1)),
       type = "l", lwd = 1.5, col = dat1Col, ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(DNAdat1, DNAranDat1, RNAdat1, RNAranDat1)))
  par(new = T)
  plot(xplot, DNAranDat2,
       ylim = c(min(DNAdat2, DNAranDat2, RNAdat2, RNAranDat2),
                max(DNAdat2, DNAranDat2, RNAdat2, RNAranDat2)),
       type = "l", lwd = 1.5, col = dat2Col, ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(DNAdat2, DNAranDat2, RNAdat2, RNAranDat2)))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = dat2Col)
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(DNAranDat1)-(flankSize/winSize), length(DNAranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(DNAranDat1)-(flankSize/winSize), length(DNAranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(DNAranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
}

# Function to plot mean coverage profile of wt vs mutant chromatin mark on same Y-axes around target and random loci
plotAvgCov_1v2v3 <- function(xplot,
                             dat1,
                             dat2,
                             dat3,
                             ranDat1,
                             ranDat2,
                             ranDat3,
                             flankSize, winSize,
                             flankLabL, flankLabR,
                             startLab1, endLab1,
                             startLab2, endLab2,
                             legendLoc, legendLabs,
                             col1, col2, col3) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3),
                max(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3)),
       type = "l", lwd = 1.5, col = col1, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, dat2, type = "l", lwd = 1.5, col = col2)
  lines(xplot, dat3, type = "l", lwd = 1.5, col = col3)
  axis(side = 2, at = pretty(c(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2, col3),
         text.col = c(col1, col2, col3),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3),
                max(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3)),
       type = "l", lwd = 1.5, col = col1, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, ranDat2, type = "l", lwd = 1.5, col = col2)
  lines(xplot, ranDat3, type = "l", lwd = 1.5, col = col3)
  axis(side = 2, at = pretty(c(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2, col3),
         text.col = c(col1, col2, col3),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot mean coverage profile of wt vs mutant chromatin mark on same Y-axes around target and random loci
plotAvgCov_1v2v3_ChIPinput <- function(xplot,
                                       dat1,
                                       dat2,
                                       dat3,
                                       ranDat1,
                                       ranDat2,
                                       ranDat3,
                                       dat4,
                                       dat5,
                                       dat6,
                                       ranDat4,
                                       ranDat5,
                                       ranDat6,
                                       flankSize, winSize,
                                       flankLabL, flankLabR,
                                       startLab1, endLab1,
                                       startLab2, endLab2,
                                       legendLoc, legendLabs,
                                       col1, col2, col3) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3,
                    dat4, dat5, dat6, ranDat4, ranDat5, ranDat6),
                max(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3,
                    dat4, dat5, dat6, ranDat4, ranDat5, ranDat6)),
       type = "l", lwd = 1.5, col = col1, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, dat2, type = "l", lwd = 1.5, col = col2)
  lines(xplot, dat3, type = "l", lwd = 1.5, col = col3)
  axis(side = 2, at = pretty(c(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3,
                               dat4, dat5, dat6, ranDat4, ranDat5, ranDat6)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2, col3),
         text.col = c(col1, col2, col3),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3,
                    dat4, dat5, dat6, ranDat4, ranDat5, ranDat6),
                max(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3,
                    dat4, dat5, dat6, ranDat4, ranDat5, ranDat6)),
       type = "l", lwd = 1.5, col = col1, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, ranDat2, type = "l", lwd = 1.5, col = col2)
  lines(xplot, ranDat3, type = "l", lwd = 1.5, col = col3)
  axis(side = 2, at = pretty(c(dat1, dat2, dat3, ranDat1, ranDat2, ranDat3,
                               dat4, dat5, dat6, ranDat4, ranDat5, ranDat6)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2, col3),
         text.col = c(col1, col2, col3),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot mean coverage profile of wt vs mutant chromatin mark on same Y-axes around target and random loci
plotAvgCov_1v2v3v4 <- function(xplot,
                               dat1,
                               dat2,
                               dat3,
                               dat4,
                               ranDat1,
                               ranDat2,
                               ranDat3,
                               ranDat4,
                               flankSize, winSize,
                               flankLabL, flankLabR,
                               startLab1, endLab1,
                               startLab2, endLab2,
                               legendLoc, legendLabs,
                               col1, col2, col3, col4) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4),
                max(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4)),
       type = "l", lwd = 1.5, col = col1, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, dat2, type = "l", lwd = 1.5, col = col2)
  lines(xplot, dat3, type = "l", lwd = 1.5, col = col3)
  lines(xplot, dat4, type = "l", lwd = 1.5, col = col4)
  axis(side = 2, at = pretty(c(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2, col3, col4),
         text.col = c(col1, col2, col3, col4),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4),
                max(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4)),
       type = "l", lwd = 1.5, col = col1, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, ranDat2, type = "l", lwd = 1.5, col = col2)
  lines(xplot, ranDat3, type = "l", lwd = 1.5, col = col3)
  lines(xplot, ranDat4, type = "l", lwd = 1.5, col = col4)
  axis(side = 2, at = pretty(c(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2, col3, col4),
         text.col = c(col1, col2, col3, col4),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot mean coverage profile of wt vs mutant chromatin mark on same Y-axes around target and random loci
plotAvgCov_1v2v3v4_ChIPinput <- function(xplot,
                                         dat1,
                                         dat2,
                                         dat3,
                                         dat4,
                                         ranDat1,
                                         ranDat2,
                                         ranDat3,
                                         ranDat4,
                                         dat5,
                                         dat6,
                                         dat7,
                                         dat8,
                                         ranDat5,
                                         ranDat6,
                                         ranDat7,
                                         ranDat8,
                                         flankSize, winSize,
                                         flankLabL, flankLabR,
                                         startLab1, endLab1,
                                         startLab2, endLab2,
                                         legendLoc, legendLabs,
                                         col1, col2, col3, col4) {
  # target loci
  plot(xplot, dat1,
       ylim = c(min(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4,
                    dat5, dat6, dat7, dat8, ranDat5, ranDat6, ranDat7, ranDat8),
                max(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4,
                    dat5, dat6, dat7, dat8, ranDat5, ranDat6, ranDat7, ranDat8)),
       type = "l", lwd = 1.5, col = col1, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, dat2, type = "l", lwd = 1.5, col = col2)
  lines(xplot, dat3, type = "l", lwd = 1.5, col = col3)
  lines(xplot, dat4, type = "l", lwd = 1.5, col = col4)
  axis(side = 2, at = pretty(c(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4,
                               dat5, dat6, dat7, dat8, ranDat5, ranDat6, ranDat7, ranDat8)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(dat1)-(flankSize/winSize), length(dat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(dat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2, col3, col4),
         text.col = c(col1, col2, col3, col4),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, ranDat1,
       ylim = c(min(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4,
                    dat5, dat6, dat7, dat8, ranDat5, ranDat6, ranDat7, ranDat8),
                max(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4,
                    dat5, dat6, dat7, dat8, ranDat5, ranDat6, ranDat7, ranDat8)),
       type = "l", lwd = 1.5, col = col1, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, ranDat2, type = "l", lwd = 1.5, col = col2)
  lines(xplot, ranDat3, type = "l", lwd = 1.5, col = col3)
  lines(xplot, ranDat4, type = "l", lwd = 1.5, col = col4)
  axis(side = 2, at = pretty(c(dat1, dat2, dat3, dat4, ranDat1, ranDat2, ranDat3, ranDat4,
                               dat5, dat6, dat7, dat8, ranDat5, ranDat6, ranDat7, ranDat8)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize), length(ranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(ranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2, col3, col4),
         text.col = c(col1, col2, col3, col4),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot mean coverage profile of wt vs mutant chromatin mark for DNA and RNA TEs on same Y-axes around target and random loci
plotAvgCov_1v2v3v4_DNA_RNA_TEs <- function(xplot,
                                           DNAdat1,
                                           DNAdat2,
                                           DNAdat3,
                                           DNAdat4,
                                           DNAranDat1,
                                           DNAranDat2,
                                           DNAranDat3,
                                           DNAranDat4,
                                           RNAdat1,
                                           RNAdat2,
                                           RNAdat3,
                                           RNAdat4,
                                           RNAranDat1,
                                           RNAranDat2,
                                           RNAranDat3,
                                           RNAranDat4,
                                           flankSize, winSize,
                                           flankLabL, flankLabR,
                                           startLab1, endLab1,
                                           startLab2, endLab2,
                                           legendLoc, legendLabs,
                                           col1, col2, col3, col4) {
  # target loci
  plot(xplot, DNAdat1,
       ylim = c(min(DNAdat1, DNAdat2, DNAdat3, DNAdat4, DNAranDat1, DNAranDat2, DNAranDat3, DNAranDat4, RNAdat1, RNAdat2, RNAdat3, RNAdat4, RNAranDat1, RNAranDat2, RNAranDat3, RNAranDat4),
                max(DNAdat1, DNAdat2, DNAdat3, DNAdat4, DNAranDat1, DNAranDat2, DNAranDat3, DNAranDat4, RNAdat1, RNAdat2, RNAdat3, RNAdat4, RNAranDat1, RNAranDat2, RNAranDat3, RNAranDat4)),
       type = "l", lwd = 1.5, col = col1, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, DNAdat2, type = "l", lwd = 1.5, col = col2)
  lines(xplot, DNAdat3, type = "l", lwd = 1.5, col = col3)
  lines(xplot, DNAdat4, type = "l", lwd = 1.5, col = col4)
  axis(side = 2, at = pretty(c(DNAdat1, DNAdat2, DNAdat3, DNAdat4, DNAranDat1, DNAranDat2, DNAranDat3, DNAranDat4, RNAdat1, RNAdat2, RNAdat3, RNAdat4, RNAranDat1, RNAranDat2, RNAranDat3, RNAranDat4)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(DNAdat1)-(flankSize/winSize), length(DNAdat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(DNAdat1)-(flankSize/winSize), length(DNAdat1)), text = c(flankLabL, startLab1, endLab1, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(DNAdat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2, col3, col4),
         text.col = c(col1, col2, col3, col4),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")

  # random loci
  plot(xplot, DNAranDat1,
       ylim = c(min(DNAdat1, DNAdat2, DNAdat3, DNAdat4, DNAranDat1, DNAranDat2, DNAranDat3, DNAranDat4, RNAdat1, RNAdat2, RNAdat3, RNAdat4, RNAranDat1, RNAranDat2, RNAranDat3, RNAranDat4),
                max(DNAdat1, DNAdat2, DNAdat3, DNAdat4, DNAranDat1, DNAranDat2, DNAranDat3, DNAranDat4, RNAdat1, RNAdat2, RNAdat3, RNAdat4, RNAranDat1, RNAranDat2, RNAranDat3, RNAranDat4)),
       type = "l", lwd = 1.5, col = col1, ann = F, xaxt = "n", yaxt = "n")
  lines(xplot, DNAranDat2, type = "l", lwd = 1.5, col = col2)
  lines(xplot, DNAranDat3, type = "l", lwd = 1.5, col = col3)
  lines(xplot, DNAranDat3, type = "l", lwd = 1.5, col = col4)
  axis(side = 2, at = pretty(c(DNAdat1, DNAdat2, DNAdat3, DNAdat4, DNAranDat1, DNAranDat2, DNAranDat3, DNAranDat4, RNAdat1, RNAdat2, RNAdat3, RNAdat4, RNAranDat1, RNAranDat2, RNAranDat3, RNAranDat4)))
  axis(side = 1, at = c(1, (flankSize/winSize)+1, length(DNAranDat1)-(flankSize/winSize), length(DNAranDat1)), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(1, (flankSize/winSize)+1, length(DNAranDat1)-(flankSize/winSize), length(DNAranDat1)), text = c(flankLabL, startLab2, endLab2, flankLabR))
  abline(v = c((flankSize/winSize)+1, length(DNAranDat1)-(flankSize/winSize)), lty = 3)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2, col3, col4),
         text.col = c(col1, col2, col3, col4),
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
