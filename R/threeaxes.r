#'  A function to plot 2 y axes
#' 
#' This function plot two response variables vs. same x axis (predictors)
#' 
#' 
#' @param data A input data frame, default to be NULL, 
#' @param eq.scale =T then y1 and y2 have same number of ticks and labels, if false, y2 scale itself.
#' @param xlab axis 1 name
#' @param y1lab axis 2 name 
#' @param y2lab axis 4 name
#' @param use is.list(x) to check if x variable is a time object...
#' @param xvar x variable, either column name or index value
#' @param yvar1 y variable, either column name or index value
#' @param yvar2 y2 variable, either column name or index value
#' @param pch1 yvar1 sympbol
#' @param col1 yvar1 color
#' @param bg1 yvar1 symbol background
#' @param pch2 yvar2 sympbol
#' @param col2 yvar2 color
#' @param bg2 yvar2 symbol background
#' @param y1log if true, yvar1 should be log transformed
#' @param y2log if true, yvar2 should be log transformed
#' @param x0lim add extra values to make sure xlim fall into the self defined range
#' @param y1lim add extra values to make sure y1lim fall into the self defined range
#' @param y2lim add extra values to make sure y2lim fall into the self defined range
#' @return Returns . 
#' @keywords respnse varialbes, three axes
#' @examples
#' xvar =1:100; yvar1 = rnorm(100)+ 1:100; yvar2 = rlnorm(100) + 1:100 
#' par(mfrow=c(1,3))
#' threeaxes.plot(my.data = NULL, xvar, yvar1, yvar2, y2log =T, eq.scale =F) 
#' threeaxes.plot(my.data = NULL, xvar, yvar1, yvar2, y2log =T, eq.scale =T) 
#' threeaxes.plot(my.data = NULL, xvar, yvar1, yvar2, y2log =F, eq.scale =F) 
#' @export
#' 
"threeaxes.plot" <- function(my.data = NULL, xvar, yvar1, yvar2, xlab = names(my.data[xvar]),
      y1lab = names(my.data[yvar1]), y2lab = names(my.data[yvar1]), pch1 = 21, col1 = 4,bg1 = "gray",
      pch2 = 22, col2 = 2, bg2 = "gray", main = "", cex.main = 1, cex.lab = 1, cex = 1, cex.axis = 1, type = "p", 
      las1 = 1, las2 = 1, xlog = FALSE, y1log = FALSE, mar = c(5,4,3,4), y2log = FALSE, eq.scale = TRUE,
      x0lim = NULL, y1lim= NULL, y2lim = NULL, ...)         {
#   par(mai =c(0,0,0,1))
   if(is.null(my.data)) {
    orig.x = xvar; orig.y1 = yvar1; orig.y2 = yvar2   
    } else {
    my.data <- my.data[order(my.data[,xvar]),]
#    good1 <- complete.cases(my.data[,xvar],my.data[,yvar1])
    orig.x <- my.data[,xvar]; orig.y1 <- my.data[,yvar1]; orig.y2 <- my.data[,yvar2]
    }
    x <- orig.x; y1 <- orig.y1; y2 <- orig.y2
    if(xlog) {
    x <-log10(orig.x); x0lim <- x0lim
    }
    if(y1log) {
    y1 <-log10(orig.y1); y1lim <- log10(y1lim); #print(y1lim); print(y1)
    }
    if(y2log) {
    y2 <-log10(orig.y2); y2lim <- log10(y2lim)
    }
    if(is.null(x0lim)) xrange = range(x)
    else xrange = range(x, x0lim)

    y1range <- range(y1, y1lim, na.rm=T)
    y2range <- range(y2, y2lim, na.rm=T)

    y2new <- (y2 - min(y2range))*diff(y1range)/diff(y2range)+ min(y1range)    # convert y2 to y1 scale

    par(mar = mar)
    plot(x,y1, xlab = xlab, ylab = y1lab, pch = pch1, col = col1, main = main, xlim= x0lim, ylim = y1lim,type=type,
        cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, bg = bg1, axes =F, ...)
    if(is.list(x)) points(as.numeric(x), y2new, col = col2, pch = pch2, type = type, bg = bg2)
    else points(x, y2new, col = col2, pch = pch2, type = type, bg = bg2)
    mtext(y2lab, side= 4, line = 2.5, cex = cex.lab) 
    box()
#### x axis
   if(xlog) {
         max.x <- ceiling(max(xrange)); min.x <- floor(min(xrange))  ## add x range for log formation
         x.at <- min.x:max.x                        # add ticks at log10 = even
         x.lab <-  round(10^x.at)                    # add labels at log10 = even number
         axis(1, at = x.at, labels = x.lab, lwd.ticks = 1.2)    # major ticks with labels
         axis(1, at = log10(1:10 * rep(x.lab[-1]/10, each = 10)), tcl = -0.3, labels = FALSE, lwd.ticks = 0.8)  # minor ticks
         x.mtick <- x.at <= max(xrange) & x.at >= min(xrange)          # major labels

         ### if only two or fewer major lables, then add tick labels at 2 and 5 log
         if(sum(x.mtick)<=2) {
           axis(1, at = log10(c(2, 5) * rep(x.lab[-1]/10, each = 2)), tcl = -0.4,lwd.ticks= 1,
           labels = (c(2, 5) * rep(x.lab[-1]/10, each = sum(x.mtick))))
          }
         }  else if(is.list(x)) axis(1, at= axTicks(1), labels = format(as.POSIXlt(axTicks(1),origin = "1970-01-01"),"%b"))
                else axis(1, las = 1)

#### y1 and y2 axis at

    if(y1log) {
         max.y1 <- ceiling(max(y1range)); min.y1 <- floor(min(y1range))  ## add y1 range for log formation
         y1.at <- min.y1:max.y1                        # add ticks at log10 = even
         y1.lab <-  round(10^y1.at)                    # add labels at log10 = even number
         axis(2, at = y1.at, labels = y1.lab, lwd.ticks = 1, col = col1, las = las1)    # major ticks with labels

         axis(2, at = log10(1:10 * rep(y1.lab[-1]/10, each = 10)) , tcl = -0.3, labels = FALSE, 
                  col = col1, lwd.ticks = 0.8)  # minor ticks
         y1.mtick <- y1.at <= max(y1range) & y1.at >= min(y1range)          # major labels
         if(sum(y1.mtick)<=2) {
           axis(2, at = log10(c(2, 5) * rep(y1.lab[-1]/10, each = 2)), tcl = -0.4,lwd.ticks= 1,
           labels = (c(2, 5) * rep(y1.lab[-1]/10, each = 2)), col = col1, las = las1)
          }
        }  else axis(2, col = col1, las = las1)
    
#########   if y2 and y1 use the same ticks 
   if(eq.scale) { 
         ### y2 axis
    
         y1new.at <- (axTicks(2) - min(y1range))*diff(y2range)/diff(y1range)+ min(y2, na.rm = T)    # convert to y2 scale
         if(y2log) {
         y2.lab <-  round(10^y1new.at)
         } else y2.lab <- round(y1new.at)    # add labels at original y2 scale
         axis(4, at = axTicks(2), labels = y2.lab, lwd.ticks = 1, col = col2, las = las2)    # major ticks with labels

   } else { ##### if y2 use different ticks as y1
      if(y2log) {
         max.y2 <- ceiling(max(y2range)); min.y2 <- floor(min(y2range))  ## add y2 range for log formation
         y2.at <- min.y2:max.y2                        # add ticks at log10 = even
         y2.lab <-  round(10^y2.at)                    # add labels at log10 = even number
         y2new.at <- (y2.at - min(y2range))*diff(y1range)/diff(y2range)+ min(y1range)    # convert to y1 scale
                  
         axis(4, at = y2new.at, labels = y2.lab, lwd.ticks = 1, col = col2, las = las2)    # major ticks with labels
         #### add minor ticks without labels
         y2.atm <- log10(1:10 * rep(y2.lab[-1]/10, each = 10))   #  desired y2 ticks
         y2new.atm <- (y2.atm - min(y2range))*diff(y1range)/diff(y2range)+ min(y1range)   # convert to y1 scale
         axis(4, at = y2new.atm, tcl = -0.3, col = col2, labels = FALSE, lwd.ticks = 0.8)  # minor ticks

        #### add pretty option if only 2 labels were added, then add more finer lables
         y2.mtick <- y2.at <= max(y2range) & y2.at >= min(y2range)          # number of major labels
         if(sum(y2.mtick)<=2) {
           y2.atm2 <- log10((c(2, 5)) * rep(y2.lab[-1]/10, each = 2))
           y2new.atm2 <- (y2.atm2 - min(y2range))*diff(y1range)/diff(y2range)+ min(y1range)
           axis(4, at = y2new.atm2, tcl = -0.4,lwd.ticks= 1, labels = (c(2, 5) * rep(y2.lab[-1]/10, each = 2)), 
           col = col2, las = las2)
         }
        }  else {
          print((y1range))
        y2.at <- (pretty(y2range) - min(y2range))*diff(y1range)/diff(y2range)+ min(y1range)
   
        axis(4, at = y2.at, labels = pretty(y2range), col = col2, las = las2)
     }   
   }
 }
