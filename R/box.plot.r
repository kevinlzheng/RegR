#' A box plot function showing pretty box plot and jitter data points
#' 
#' This function box plot similar to ggplot2 function
#' modified from https://github.com/zonination/perceptions/blob/master/percept.R)
#' 
#' @param x A data frame or matrix showing with cross tabled format, 
#' @param horizontal a boolean variable deciding if the box should be horizontal 
#' @param outline if an outline 
#' @param mar margin for plot area
#' @param xlab
#' @param ylab
#' @param main for title
#' @param names group names, can be entered manually
#' @param line when mar changes, ylab should be changed to by line = 3
#' @examples
#' set.seed(100)
#' n <- rnorm(11, 2, 11)
#' x <-  (matrix(rnorm(15*11, 3, 7), nrow=15, ncol=11) *2)
#' x2 <- sapply(1:11, function(i) x[,i]+n[i]) 
#' names <- paste("Var", 1:11, sep =".")
#' box.plot (x2, horizontal = TRUE, outline = FALSE, lty=1, staplewex=0, 
#'          boxwex=0.8, boxlwd=1, medlwd=1,xaxt="n", yaxt="n", notch = FALSE,
#'          mar=c(5,5,3,1), xlab ="Assigned Probability %", ylab = "Phrase",
#'          main = "Perceptions of Probability", names = names) 
#' @export

box.plot <- function(x, horizontal = TRUE, outline = FALSE, lty=1, staplewex=0, line=6, 
                     boxwex=0.8, boxlwd=1, medlwd=1,xaxt="n", yaxt="n", notch = FALSE,
                     mar=c(5,8,3,1), xlab ="Assigned Probability %", ylab = "Phrase",
                     main = "Perceptions of Probability",cex.main = 1,names = NULL) {
   # load the data and convert to a matrix
    m <- as.matrix(x[,ncol(x):1])
    yrange <- range(m, na.rm =T)
    ypretty <- pretty(yrange)
    npretty <- length(ypretty)
    ymin = ypretty[1]
    ymax = ypretty[npretty]
    ydiff <- ymax - ymin
   # create some random data for jitter
    r <-  (matrix(runif(nrow(m)*ncol(m)), nrow=nrow(m), ncol=ncol(m)) / 2) - 0.25

   # create colours and colour matrix (for points)
    cols  <- colorRampPalette(brewer.pal(12, "Set3"), alpha=TRUE)(ncol(m))
    colsm <- matrix(rep(cols, each=nrow(m)), ncol=ncol(m))

   # get the greys
    palette <- brewer.pal("Greys", n=9)

   # set graphical area
    par(bty="n", bg=palette[2], mar = mar)
   
   # plot initial boxplot
    boxplot(m~col(m),horizontal=horizontal, outline=outline, lty=lty, staplewex=staplewex, ylab = ylab,
            boxwex=boxwex, boxlwd=boxlwd, medlwd=medlwd, col=cols, xaxt=xaxt, yaxt=yaxt, add =F)
#    grid(nx = npretty, ny = ncol(m), col = palette[5], lty = 1)
   # plot gridlines
    for (i in seq(ymin, ymax,length = npretty)) {
       lines(c(i,i), c(0,ncol(m)+1), col=palette[5])
       }

    for (i in seq(1,ncol(m),by=1)) {
       lines(c(ypretty[1]-0.05*ydiff, ypretty[npretty]+0.05*ydiff), c(i,i), col=palette[5])
       }

    # plot points
    points(m, col(m)+r, col=colsm, pch=16)
    ynames <- names
    if(is.null(names)) ynames <- colnames(m)

    # overlay boxplot
    boxplot(m~col(m), horizontal=horizontal, outline=outline, lty=lty, staplewex=staplewex, 
            boxwex=boxwex, boxlwd=boxlwd, medlwd=medlwd, col=cols, xaxt=xaxt, yaxt=yaxt, add =T)

    # add axes and title
    axis(side=1, at= ypretty, col.axis=palette[7], cex.axis=0.8, lty=0, tick=NA, line=-1)
    axis(side=1, at= mean(ypretty), labels= xlab, lty=0, tick=NA, col.axis=palette[7])
    axis(side=2, at=1:ncol(m), col.axis=palette[7], cex.axis=0.8, lty=0, tick=NA, labels= ynames, las=2)
    axis(side=2, at= ncol(m)/2, labels=ylab, col.axis=palette[7], lty=0, tick=NA, las=3, line= line)
    title(main, cex = cex.main)
}


