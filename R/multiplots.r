#' plot multiple color graphs
#' 
#' The output of this function to return 1.Weighted Average, 2. cdf_Abundance based, 3. cdf_ presence/absence based;
#' 4. ecdf weighted, 5. cdf weight new; 6. Linear logistic regression, 7. quadratic logistic 8. GAM  5~7 using full data range; 
#' 9~11. repeat 6~8 but uses observed range for each single taxon; 12 Count. 13. Raw quantiles
#'
# KevinLZheng2000@gmail.com 
#################
#' @param my.data data frame 
#' @param xvar x variable, either column name or index value
#' @param yvar y variable, either column name or index value
#' @param fit fitted line, either loess or linear or none
#' @param classes group and a factor variable
#' @param interval "prediction" or "confidence" interval
#' @param cols colors for different classes
#' @param ltys line types for different groups
#' @param pchs point type 
#' @param add.ln add fitted lines
#' @param predx predicted value at 100%
#' @param add.pt add predicted point
#' @keywords scatter plot, interval, change point, regression
#' @examples
#' set.seed(1)
#' n = 50
#' x <- rlnorm(n) + 1:n
#' x <- (x - min(x))/max(x)
#' y <- n + 1:n + (rnorm(n))*0.2*n
#' g <- paste0("g", sample(1:5, n))
#' par(mfrow=c(1,2))
#' plot(y~x)
#' multiplots(data = data.frame(x,y, g), xvar= "x", yvar = "y", classes = "g", fit = "linear")
#' @export

multiplots <- function (my.data = my.data, xvar, yvar, classes, xlab = names(my.data[xvar]),
    fit = "loess", log = "", ylab = names(my.data[yvar]), type = "n", pos = "topleft", 
    cols = 1:20, pchs = 1:20, rounder1 = 2, rounder2 = 2, ltys = 1:20, main = "", 
    rounder = 2, interval = "prediction", level = 0.5, predx = 100, add.pt = TRUE,
    add.ln = FALSE, lty = 1, cex = 1, las = 1, cex1 = 1,...) {

    good <- complete.cases(my.data[,xvar],my.data[,yvar],my.data[,classes])
    dataset <- my.data[good,]; print(dim(dataset))
    my.class <- levels(dataset[,classes])
    len <- length(my.class)
    x.all <- dataset[,xvar]; y.all <- dataset[,yvar]
      if(log =="x" | log == "xy") {
          x.all <- log10(dataset[,xvar])
          if(predx !="") predx = log10(predx)
          }
      if(log =="y" | log == "xy") {
          y.all <- log10(dataset[,yvar])  
          } 
          
    plot(x.all, y.all, type = type, xlab = xlab, ylab = ylab, axes = F, main = main, ...)
    newset <- data.frame(x = x.all, y = y.all, classes = dataset[,classes]) ; print(len)
    for (i in 1:len) {
     my.subset <- subset(newset, classes == my.class[i]) ;print(dim(my.subset))
          print(my.class[i])
         newx <- seq(min(my.subset$x), max(my.subset$x), length =100)
        if(add.pt) {
         points(my.subset$x, my.subset$y, col = cols[i], pch = pchs[i], cex = cex1) 
        } 
        if(add.ln) {
         lines(my.subset[order(my.subset$x),"x"], my.subset[order(my.subset$x),"y"], col = cols[i], lty = ltys[i]) 
        }
        #### add model
        if(nrow(my.subset)> 20) {
            if(fit == "loess") {
               model <- loess(y~x, data = my.subset)
               pred <- predict(model, newdata = data.frame(x=newx), se.fit=T)
               points(newx, pred, pch = pchs[i], col = cols[i], type = "l", lty = ltys[i])
               }
            ##### linear model add model lines
            if(fit == "linear") {   
               model <- lm(y~x, data = my.subset)
               pred <- predict(model, newdata = data.frame(x=newx), interval = interval, level = level, se.fit=T)
                #### add prediction lines
               if(interval!="none") {
                points(newx, pred$fit[,1], pch = pchs[i], col = cols[i], type = "l", lty = 1)
                points(newx, pred$fit[,2], pch = pchs[i], col = cols[i], type = "l", lty = 2)
                points(newx, pred$fit[,3], pch = pchs[i], col = cols[i], type = "l", lty = 2)
                }  else {
                points(newx, pred$fit, pch = pchs[i], col = cols[i], type = "l", lty = lty[i]); print(cols[i])
                }
              ## interpolation of max or min prediction
              if(predx!="") {
               predy <- predict(model, newdata = data.frame(x=predx), interval = interval, level = level, se.fit=T)    
               print(predy)
               ### alternative calcualtion
               if(interval!="none") {
               pred.up.pr <- predy$fit[,1] + qt((1-level)/2,13)*sqrt(predy$se.fit^2+ predy$resid^2)
               pred.lo.pr <- predy$fit[,1] + qt((1+level)/2,13)*sqrt(predy$se.fit^2+ predy$resid^2)
               points(rep(predx,3), c(predy$fit[,1], pred.up.pr, pred.lo.pr), col = cols[i], pch = 19)
               cat("Region", my.class[i],c(predx, predy$fit[,1], pred.up.pr, pred.lo.pr), "\n")
               if(log =="y" | log == "xy") {
                 mtext(round(exp(predy$fit[,1]), rounder), side = 4, at = c(100, predy$fit[,1]), las = 1, col = cols[i])
                 cat("Region", my.class[i],predx, round(exp(c(predy$fit[,1], pred.up.pr, pred.lo.pr)), rounder), "\n")
                } else {
                 mtext(round(predy$fit[,1], rounder), side = 4, at = c(100, predy$fit[,1]), las = 1,col = cols[i])
                 cat("Region", my.class[i],c(predx, predy$fit[,1], pred.up.pr, pred.lo.pr), "\n")                 
                }
               }
             }  # end interpolation
           } ####end of linear model

      }  ##### end add  model 

    }
    if(log == "x" |log == "xy") {
         max.pow <- ceiling(max(newset$x)); min.pow <- floor(min(newset$x))  ## add x range for log formation
         axis.at <- min.pow:max.pow                        # add ticks at log10 = even
         axis.lab <-  round(10^axis.at, rounder1)                    # add labels at log10 = even number
         axis(1, at = axis.at, labels = axis.lab, lwd.ticks = 1, las = las)    # major ticks with labels
         axis(1, at = log10(1:10 * rep(axis.lab[-1]/10, each = 10)), tcl = -0.3, labels = FALSE, lwd.ticks = 0.8)  # minor ticks
         xtick <- axis.at <= max(newset$x) & axis.at >= min(newset$x)          # major labels
        
         ### if only two or fewer major lables, then add tick labels at 2 and 5 log
         if(sum(xtick)<=2) {
           axis(1, at = log10(c(2, 5) * rep(axis.lab[-1]/10, each = sum(xtick))), tcl = -0.4,lwd.ticks= 0.8, 
           labels = (c(2, 5) * rep(axis.lab[-1]/10, each = sum(xtick))), las = las)
          }
       }  else axis(1,las = las)
        
    if(log == "y" |log == "xy") {
         max.pow <- ceiling(max(newset$y)); min.pow <- floor(min(newset$y))  ## add x range for log formation
         axis.at <- min.pow:max.pow                        # add ticks at log10 = even
         axis.lab <-  round(10^axis.at,rounder2)                    # add labels at log10 = even number
         axis(2, at = axis.at, labels = axis.lab, lwd.ticks = 1, las = las)    # major ticks with labels
         axis(2, at = log10(1:10 * rep(axis.lab[-1]/10, each = 10)), tcl = -0.3, labels = FALSE, lwd.ticks = 0.8)  # minor ticks
         xtick <- axis.at <= max(newset$y) & axis.at >= min(newset$y)          # major labels
        
         ### if only two or fewer major lables, then add tick labels at 2 and 5 log
         if(sum(xtick)<=2) {
           axis(2, at = log10(c(2, 5) * rep(axis.lab[-1]/10, each = sum(xtick))), tcl = -0.4,lwd.ticks= 0.8, 
           labels = (c(2, 5) * rep(axis.lab[-1]/10, each = sum(xtick))), las = las)
          }
       }  else axis(2,las = las)
    if(pos!="") {    
    legend(pos, pch = pchs[1:len], col = cols[1:len], lty = lty, legend= c(my.class), title = classes, cex = cex,...)
    }
    box(bty="l")
  }
