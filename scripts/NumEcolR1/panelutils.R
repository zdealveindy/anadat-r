# panelutils.R 
#
# License: GPL-2 
# Author: Francois Gillet, February 2007
#
  ## Put Pearson, Spearman or Kendall correlations on the upper panel
     panel.cor <- function(x, y, method="pearson", digits=3, cex.cor=1.2)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- cor(x, y, method=method)
         ra <- cor.test(x, y, method=method)$p.value
         txt <- round(r, digits)
         sig <- 1
         prefix <- ""
         if(ra <= 0.1) prefix <- "."
         if(ra <= 0.05) prefix <- "*"
         if(ra <= 0.01) prefix <- "**"
         if(ra <= 0.001) prefix <- "***"
         if(ra <= 0.001) sig <- 2
         color <- 2
         if(r < 0) color <- 4
#         color <- "gray10"
#         if(r < 0) color <- "gray50"
         txt <- paste(txt, prefix, sep="\n")
         text(0.5, 0.5, txt, cex = cex.cor, font=sig, col=color)
     }

  ## Put histograms on the diagonal
     panel.hist <- function(x, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(usr[1:2], 0, 1.5) )
         h <- hist(x, plot = FALSE)
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
#         rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
     }


#Usage:
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, method="kendall")
