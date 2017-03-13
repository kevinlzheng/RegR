#' A regression function for simple regression and break point analysis and plotting
#' 
#' This function performed basic regression analysis with the option to add confidence intervals and
#' prediction intervals; alternative regression approaches, e.g., segmented regression, change point
#' analysis (through rpart) and confidence bounds, quantile regression, and lowess regression splines.
#' This version also added a capacity for transforming percent values (where 0 is probably the minimum)
#' by adding 1 or 100 before log transformation if users choose so. 
#' Any x variable with 0s and log is true, the add 1 to the original or 100 times of the original.
#' 
#' @param data A input data frame, default to be NULL, 
#' @param x.add 1 to convert percent x value = 0 to 1+ x value , or 100, then *100+1, 
#' @param xName xName for making equation labels
#' @param yName yName for making equation labels
#' @param sign if TRUE then use signif for equation otherwise, use round
#' @param xvar x variable, either column name or index value
#' @param yvar y variable, either column name or index value
#' @param yvar for column index of x, y values in data data frame, or column names
#' @param add.fit either "linear", "lowess", "loess", "segment" (segmented regression), "change" (regression tree)
#' @param inteval "confidence", "prediction", "none" 
#' @param add.cor add spearman correlation, 
#' @param add.reg will add regression model
#' @param log will control x,y axis if log, 
#' @param rq choose if we want to do quantile regression
#' @param rq.tau a vector to indicate which quantiles we want to use
#' @param tx.col is for col of labs ln.col for fitted line color
#' @param lab.pos control position of labels "bottomright" "bottomleft", "topright", "topleft"
#' @param add.r2 a boolean to indicate if r2
#' @param nlrq only works for logistic fit
#' @param add.interp will add a interpolation for inverse prediction of x based on pred.y
#' @param add if TRUE, then add to existing graph
#' @param addaxes to plot axes
#' @param bty add a box frame
#' @param xlog.scale either 10 or exp
#' @param ylog.scale either 10 or exp
#' @param xlim xlim
#' @param ylim ylim
#' @param xlab xlabel
#' @param ylab ylabel
#' @param addaxes if true, add axes
#' @param ex.eq equation
#' @param cex  magnifier
#' @param cex.axis # magnifier
#' @param pch point type
#' @param lty line type
#' @param las direct
#' @param bg background col
#' @param side which side to add
#' @param x.add add x
#' @param labs add labels
#' @param add.rq # add regression quantile
#' @param rounder1 round number of digits for x value
#' @param rounder2 round number of digits for y value
#' @param pred.x used to predict x value
#' @param pred.y used to predict y 
#' @param level CI level
#' @return Returns . 
#'  sum((x-xmean)^2)) = sum(x^2)- (sum(x))^2/length(x)
#   regression predicted SE = residual error * sqrt(1/n + (x0-xmean)^2/sum((x-xmean)^2)) 
##  regression prediction for a single value SE = residual error * sqrt(1/n + (x0-xmean)^2/sum((x-xmean)^2)) 
## residual SE example   sqrt(((x[1]-mean(x))^2/sum((x-mean(x))^2)+1/length(x)))*pred.line$residual.scale
#' @keywords scatter plot, interval, change point, regression
#' @examples
#' set.seed(1)
#' n = 50
#' x <- rlnorm(n) + 1:n
#' x <- (x - min(x))/max(x)
#' y <- n + 1:n + (rnorm(n))*0.2*n
#' par(mfrow=c(1,2))
#' plot(y~x)
#' scatter.plot(data = data.frame(x,y), xvar= "x", yvar = "y", add.fit = "linear", xlim = c(0,n), 
#'     ylim = range(y), log = "x", x.add =100)
#' @export
#' 
#' 

"scatter.plot" <- function(data = NULL, xvar, yvar, xName = "x", yName = "y", 
      xlim = NULL, ylim = NULL, sign = FALSE,
      xlab = names(data[xvar]), ylab = names(data[yvar]),  
      addaxes= TRUE, cex.eq = 1.5, bty="o", type = "p", pch = 1, col = 1, 
      main = "", cex.main = 1, cex.lab = 1, cex = 1, cex.axis = 1,  lty = 1, 
      las = 1, bg= "gray", side = 4, tx.col=2,  xlog.scale = 10, ylog.scale =10, 
      lab.pos = "topleft", labs = "", x.add = 0, 
      add.fit = "none", interval = "none", ln.col = "red", level = 0.95, 
      add.cor = FALSE, add.reg = FALSE, add.rq = FALSE, 
      nl.rq = FALSE, rq.tau = 0.5, log = "",  add = FALSE, add.r2 = FALSE,
      add.interp = FALSE, pred.x = NA, pred.y = NA,  
      rounder1 = 2, rounder2 = 0,...)
{
######################## preparing data    
    if(is.null(data)) {
    orig.x <- xvar[!is.na(xvar) & !is.na(yvar)]
    orig.y <- yvar[!is.na(xvar) & !is.na(yvar)]
    x = orig.x; y = orig.y
    } else {
    data <- data[order(data[,xvar]),]
    good <- complete.cases(data[,xvar],data[,yvar])
#     options(scipen=3) 
    orig.x <- data[,xvar][good]; orig.y <- data[,yvar][good]
    x <- orig.x ; y<- orig.y
    }

    if(is.null(xlim))  xlim <- range(x)
    else             xlim <- xlim
    if(is.null(ylim) ) ylim <- range(y)  
    else             ylim <- ylim    

    # log transform values
    if(log == "x" |log == "xy") {
       if(x.add == 0) {
          x <-log(x, base = xlog.scale)
          xlim <- log(xlim, base = xlog.scale)
         } else {                    # here is where you have different log(x) value, add 1 
           if (x.add ==1) {
            x <-log(x + 1, base = xlog.scale)
           } 
           if(x.add == 100){
            x <-log(x*100 + 1, base = xlog.scale)          
           }
            xlim <- log(xlim + 1, base = xlog.scale)  
         }
    }
    if(log =="y" |log == "xy") {
          y <- log(y,  base = xlog.scale)
          ylim <- log(ylim, base = xlog.scale)
     }
    ## Set up dummy x values for plotting
    xrange = range(x)
    yrange = range(y)

    xlims <- range(xlim, xrange) 
    ylims <- range(ylim, yrange)

    newx = seq(xrange[1], xrange[2], length=1000)
    newdata = data.frame(x = newx)
############### legend function
  
cor.pos <- function(x, y,...) {
    xtext <- c("bottomright","bottomleft", "topleft", "topright")
    box1 <- rep(NA, length(xtext))
    for(i in 1:length(xtext)) {
         xx <- legend(xtext[i],plot=FALSE,legend ="", ...)
         x1 <- xx$rect$left; x2 <- xx$rect$left+xx$rect$w
         y2 <- xx$rect$top; y1 <- y2-xx$rect$h
         box1[i] <- sum(x>= x1 & x <=x2 & y >= y1 & y <= y2)
      }
   return(xtext[which.min(box1)][1])
  }
    if(lab.pos =="") lab.pos <- cor.pos(x,y)
#### Draw the plot with appropiate axes 
    if(!add) {
    plot(x, y, xlim = xlims, ylim = ylims, axes = FALSE,main = main, 
         xlab = xlab, ylab = ylab, type = type, pch = pch, col = col,
         bg = bg, cex.main = cex.main, cex.lab = cex.lab, cex = cex, 
         cex.axis = cex.axis, ...) 
     } else {
     points(x, y, pch = pch, type = type, col = col, bg = bg, cex = cex)
     }    
    xmin <- ifelse(min(xlims) < log10(0.20*10^ceiling(min(xlims))), 
                   floor(min(xlims)), min(xlims)) 
    ymin <- ifelse(min(ylims) < log10(0.20*10^ceiling(min(ylims))), 
                   floor(min(ylims)), min(ylims))
      
############################ Calculate linear regression, and obtain range for y-axis
    if(add.fit != "linear" & add.cor) {
      correl <- cor(x,y, method = "spearman")
      legend(lab.pos, bty = "n", legend = bquote(italic(r) ==.(round(correl, digit = 2))))
      } 
    
    if(add.fit == "linear"){
      my.fit = lm(y ~ x)    
      pred.line = predict(my.fit,newdata)             
      yrange = range(c(yrange,pred.line),na.rm=TRUE)
      if(interval =="confidence"){
        conf.fit = predict(my.fit,newdata,interval="confidence",level = level)
        yrange = range(c(yrange,conf.fit[,2:3]),na.rm=TRUE)
      }
      if(interval =="prediction"){
        pred.fit = predict(my.fit,newdata,interval="prediction",level = level)
        yrange = range(c(yrange,pred.fit[,2:3]),na.rm=TRUE)
        if(!is.na(pred.x)) {
          pred2 <- predict(my.fit, newdata = data.frame(x = pred.x), 
                   interval="prediction",level= level)
          print(10^pred2)
        }
        if(add.interp) {
          print(anova(my.fit))
        ss = anova(my.fit)["Residuals","Sum Sq"] 
        df = anova(my.fit)["Residuals","Df"]
        intcpt <- (pred.y - my.fit$coef[1])/my.fit$coef[2]

        hi.int <- intcpt + sqrt(ss/df)/my.fit$coef[2] * qt((1+ level)/2, df)
        lo.int <- intcpt + sqrt(ss/df)/my.fit$coef[2] * qt((1- level)/2, df)

        #### alternative methods for prediction

        out.search <- c(newx[which.min(abs(pred.fit[,1]-pred.y))],
                        newx[which.min(abs(pred.fit[,2]-pred.y))],
                        newx[which.min(abs(pred.fit[,3]-pred.y))])    
                  if(log =="x"|log=="xy") {
                    out.search = 10^out.search
                  }
        print(cat("interpolation for inverse prediction of x based on pred.y\n",out.search,                  
              "This is search results")) 
        }
      }
    }   

    ## Calculate loess smooth if requested, and obtain range for y-axis

    if(add.fit == "loess"){
      lo.fit = loess(y~x)
      
      lo.pred.line = predict(lo.fit,newdata,se=TRUE)
      yrange = range(c(yrange,lo.pred.line$fit),na.rm=TRUE)
      if(interval == "confidence"){
        lo.conf.shift = lo.pred.line$se.fit * qt( (1+level)/2, lo.pred.line$df )
        lo.conf.low = lo.pred.line$fit - lo.conf.shift
        lo.conf.high = lo.pred.line$fit + lo.conf.shift
        yrange = range(c(yrange,lo.conf.low,lo.conf.high),na.rm=TRUE)
      }
      if(interval == "prediction"){
        lo.pred.shift = sqrt(lo.pred.line$se.fit^2 + lo.pred.line$residual.scale^2 ) *
          qt( (1+level)/2, lo.pred.line$df )
        lo.pred.low = lo.pred.line$fit - lo.pred.shift
        lo.pred.high = lo.pred.line$fit + lo.pred.shift
        yrange = range(c(yrange,lo.pred.low,lo.pred.high),na.rm=TRUE)
      } 
    }                          

    if(add.fit == "segment"){  # my own function
      x.sort <- sort(unique(x))
      tmp.res <- rep(NA, length(x.sort)-6)
      for(i in 4:(length(x.sort)-3)) {    # find all possible segments
      x1 <- x[x<=x.sort[i]]; y1<- y[x<=x.sort[i]]
      x2 <- x[x>x.sort[i]]; y2 <- y[x>x.sort[i]]
      mod1 <- lm(y1~x1); mod2<- lm(y2~x2)
      tmp.res[i-3] <- anova(mod1)[2,3] + anova(mod2)[2,3]    # total residual
      }
      seg.pt <-   x.sort[which.min(tmp.res)+3]         # minimum residual
      if(x.add!=1 & x.add!=100) {
      iprint <-ifelse(log =="x"|log =="xy", round(xlog.scale^seg.pt, rounder1),
                      round(seg.pt,rounder1))
      } else {
      iprint <-ifelse(log =="x"|log =="xy", round(xlog.scale^seg.pt-1, rounder1),
                      round(seg.pt,rounder1))        
      }
      print(paste("segment",iprint))  
      Ind <- ifelse(x <= seg.pt, 1, 0)
      seg.fit = lm(y~x + I((x-seg.pt)*Ind))
      pred.y <- predict(seg.fit, newdata = data.frame(x = x.sort, 
                        Ind = x.sort<=seg.pt ) )
      mtext(iprint, line = 0.2, side =3)
      lines(x.sort[x.sort<=seg.pt], pred.y[x.sort<=seg.pt], lty = 2, col = ln.col, lwd = 1.5)
      lines(x.sort[x.sort>seg.pt], pred.y[x.sort>seg.pt], lty = 2, col = ln.col, lwd = 1.5)
    }
    if(add.fit == "segmented") {
      library(segmented)
      linear.model <- lm(y~x)
      seg.model <- segmented(linear.model,seg.Z=~x,psi = median(x),#  quantile(x, prob = c(0.25, 0.75)), 
                             control=seg.control(display=FALSE))
   #   if(nrow(confint(seg.model)[[1]])==1 {
      brk.pt <- confint(seg.model)[[1]][1:3]      
      if(x.add!=1 & x.add!=100) {
      iprint <-ifelse(log =="x"|log =="xy", round(xlog.scale^brk.pt, rounder1),
                      round(brk.pt,rounder1))
      } else {
      iprint <-ifelse(log =="x"|log =="xy", round(xlog.scale^brk.pt-1, rounder1),
                      round(brk.pt,rounder1))        
      }
      rect(brk.pt[2], min(y),brk.pt[3], max(y), col = "#11111120", border =NA)
    #  print(summary(seg.model)) 
    #  print(brk.pt)
      }
   if(!add) {
   if(addaxes) {           # add axes 
    # for axis 1 mostly transformed data
    if(log == "x" |log == "xy") {
       if(x.add==0) {
         max.pow <- ceiling(log10(xlog.scale^(max(xlims))))
         min.pow <- floor(log10(xlog.scale^(min(xlims))))          
         axis.lab <- round(10^(min.pow:max.pow), rounder1)
         axis.at <- log(axis.lab, xlog.scale)    # add ticks at log10 = even
         } else {
           max.pow <- ceiling(max(xlims))
           min.pow <- floor(min(x[x != 0]))        # define min landuse 0.01% (log10(0.0001*100+1))=           
           axis.lab <- round(10^(min.pow:max.pow), rounder1)
           axis.at <- log(axis.lab +1, xlog.scale) # add ticks at log10 = even 
         }
         axis(1, at = axis.at, labels = axis.lab, tcl = -.5, las = las ) 
            xtick <- axis.at <= max(xlims) & axis.at >= min(xlims)          # major labels
  ### if only two or fewer major lables, then add tick labels at 2 and 5 log
         mark1 <- 5; mark2 <- c(2,5)
              if(x.add == 1 | x.add ==100) { mark1 <- 6; mark2 <- c(3, 6) }
      if(sum(xtick)==2) {
              at1 <- log(mark1 * rep(axis.lab[-1]/10, each = 1), base =  xlog.scale)
              axis(1, at = at1, tcl = -0.4, lwd.ticks= 1, labels = (5 * rep(axis.lab[-1]/10, each = 1)), las = las)
          }   else if (sum(xtick)<2) {
              at2 <- log(mark2 * rep(axis.lab[-1]/10, each = 2),  base = xlog.scale)
              axis(1, at = at2, tcl = -0.4, lwd.ticks= 1, labels = (c(2, 5) * rep(axis.lab[-1]/10, each = 2)), las = las)
          }            
          if(sum(xtick) < 6) {
              at3 <- log(1:10 * rep(axis.lab[-1]/10, each = 10), base = xlog.scale)    
              axis(1, at = at3, tcl = -0.3, labels = FALSE, lwd.ticks = 0.8)
          }  
    }  else axis(1,las = las)     ## nontransformed data
    ##### for axis 2 transformed data 
    if(log=="y" |log=="xy") {
         max.pow <- ceiling(log10(ylog.scale^(max(ylims))))
         min.pow <- floor(log10( ylog.scale^(min(ylims))))
       
         axis.lab <- round(10^(min.pow:max.pow), rounder2)
         axis.at <- log(axis.lab, ylog.scale) 
         ytick <- axis.at <= max(ylims) & axis.at >= min(ylims)

         axis(2, at = axis.at, labels = axis.lab, tcl = -.5, las = las)
         if(sum(ytick)== 2) {
            at1 <- log(5 * rep(axis.lab[-1]/10, each = 1), base = ylog.scale)
            axis(2, at = at1, tcl = -0.4, lwd.ticks= 1, labels = ( 5 * rep(axis.lab[-1]/10, each = 1)), las = las)
           }   else if (sum(ytick) < 2) { 
            at2 <- log(c(2, 5) * rep(axis.lab[-1]/10, each = 2), base = ylog.scale)          
            axis(2, at = at2, tcl = -0.4, lwd.ticks= 1, labels = (c(2, 5) * rep(axis.lab[-1]/10, each = 2)), las = las)    
           }  
           if(sum(ytick) <6) {
              at3 <- log(1:10 * rep(axis.lab[-1]/10, each = 10), base = ylog.scale)    
              axis(2, at = at3, tcl = -0.3, labels = FALSE, lwd.ticks = 0.8)
          }  
      } else axis(2, las = las)
      box(bty = bty) 
     } 
    }
    if(add.fit == "linear"){
      lines(newx,pred.line, lty = lty, col = ln.col, lwd = 1.5)
      
      ## Add confidence interval lines
      if(interval == "confidence"){
        lines(newx,conf.fit[,2],lty = 2,col = ln.col, lwd = 1.5)
        lines(newx,conf.fit[,3],lty = 2,col = ln.col, lwd = 1.5)
      }
      
      ## Add prediction interval lines
      if(interval == "prediction"){
        lines(newx,pred.fit[,2],lty = 2,col = ln.col)
        lines(newx,pred.fit[,3],lty = 2,col = ln.col)
        if(add.interp) {
          points(c(intcpt, lo.int, hi.int), rep(pred.y, 3), col = "red", cex = 2)
          if(log =="x"|log == "xy") {
            out.int = c(10^intcpt, 10^lo.int, 10^hi.int)
          }
          print(cat("interpolation for inverse prediction of x based on pred.y:\n", out.int))
        }
      }

      ## Add regression info
      my.oper = ifelse(my.fit$coef[1] > 0, " + ", " - ")
      equation <- paste(yName, " = ", round(my.fit$coef[2],3), "*", xName,
                 my.oper, round(abs(my.fit$coef[1]),3), sep="")
      if(sign) {
        equation <- paste(yName, " = ", signif(my.fit$coef[2],3), "*", xName,
                          my.oper, signif(abs(my.fit$coef[1]),3), sep="")
      }
      rsq <- bquote(italic("R"^{2}) ==.(signif(summary(my.fit)$r.squared,2)) )

      if(add.r2) {
        legend(lab.pos, bty = "n", legend = rsq)
      } 
      # else mtext(rsq, side=side, line=0.2)
      if(add.reg) {
        if(my.fit$coef[2]>0) {
          if(lab.pos == "topleft") {
          text(max(x), min(y), equation, pos=2, cex = cex.eq)
        } else {
          text(min(x), max(y), equation, pos=4, cex = cex.eq)
        } 
      }

      # add regression equation in the plot
      if(my.fit$coef[2]<0) {
        if(lab.pos == "topright") {
        text(min(x), min(y), equation, pos = 4, cex = cex.eq)
        } else {
        text(max(x), max(y), equation, pos = 2, cex = cex.eq)
        }
      }
     }
    }

      ## Add Change point analysis
    if(add.fit == "change")  { 
      pvalue <- chngp.nonpar(data, xvar, yvar,nindex=NULL, num = T)

#      if(interval == "confidence") {
        theta <- function(i) {rpart(y[i] ~ x[i])$splits[1, 4]}
        if (!is.null(theta)) {
	       results <- bootstrap(1:length(x), 500,theta)
	       bpoint <- unlist(results$thetastar)
         ci <- quantile(bpoint, p=c((1 - level)/2, 0.5,(1 + level)/2)) 
         }
  
        rect(ci[1], min(y),ci[3], max(y), col = "#11111120", border =NA)

		    if(log == "x" | log == "xy") {
          if(x.add!=1 & x.add!=100) iprint <- round(xlog.scale^ci, rounder1)
          else iprint <- round(xlog.scale^ci - 1 , rounder1)

         } else iprint <- round(ci, rounder1)
      
   	      print(paste(iprint[1], iprint[2], iprint[3], sep ="~"))  
		    #  abline(v = ci[2], lty = lty, lwd = 1.5, col = ln.col)
#        }       
#      return(pvalue)
      }
    ## add lowess lines
    if(add.fit == "lowess") {
        lines(lowess(x, y, f = 0.75), lty = lty, col = ln.col)
    }  

    ## Draw loess output
    if(add.fit == "loess" ){
      lines(newx,lo.pred.line$fit,lty = lty,col = ln.col)
      
      ## Add confidence interval lines
      if(interval == "confidence"){
        lines(newx, lo.conf.low, lty = 2, col = ln.col)
        lines(newx, lo.conf.high, lty = 2,col = ln.col)
      }
      
      ## Add prediction interval lines
      if(interval == "prediction"){
        lines(newx, lo.pred.low, lty=3, col= ln.col)
        lines(newx, lo.pred.high, lty=3, col= ln.col)
      }
    }

      # fit a quantile regression model
     if(add.rq) {
         for(i in 1:length(rq.tau)) {
          if(nl.rq) { 
          my.rq <- nlrq(as.formula(y ~ SSlogis(x, Asym, mid, scal)), trace=FALSE, tau = rq.tau[i])        
          print(my.rq)
          } else  {          
          my.rq <- rq(as.formula(y ~ x), tau = rq.tau[i])
          }
          yy <- predict(my.rq, newdata = newdata, interval = "confidence", level = level)
          yrange = range(c(y,yy,na.rm=TRUE))
          cols <- ifelse(length(rq.tau)==1, ln.col, 16 + i)
          if(nl.rq) {
          lines(newdata$x, yy, col = cols, lty = "solid", lwd = 1.5 )
          } else {
          lines(newdata$x, yy[,1], col = cols, lty = "solid", lwd = 1.5 )
          }
      ## Add regression info
      if(add.reg) {
          my.oper = ifelse(my.rq$coef[1] > 0, " + ", " - ")
          equation <- paste(yName, " = ", signif(my.rq$coef[2],3), "*", xName,
                 my.oper, signif(abs(my.rq$coef[1]),3), sep="")          
        if(my.rq$coef[2]>0) {
          if(lab.pos == "topleft") {
          text(max(x), min(y), equation, pos=2, cex = cex.eq)
          } else {
          text(min(x), max(y), equation, pos=4, cex = cex.eq)
           } 
         }
        if(my.rq$coef[2] < 0) {
            if(lab.pos == "topright") {
            text(min(x), min(y), equation, pos = 4, cex = cex.eq)
            } else {
            text(max(x), max(y), equation, pos = 2, cex = cex.eq)
           }
         }          
        }    ## output coefficient CIs if necessary
          if (!(nl.rq) & interval !="none") {
          print(paste("Quantile Confidence Interval Tables for quantile", rq.tau[i]))
          print(summary(my.rq))
          yrange = range(c(yrange,yy[,2],yy[,3]),na.rm=TRUE)
          lines(newdata$x, yy[,2], col = 16 + i, lty = "dashed", lwd = 1.5 )
          lines(newdata$x, yy[,3], col = 16 + i, lty = "dashed", lwd = 1.5 )
          } else {
           ## output just the coefficient estimates
          print(paste("Quantile Coefficient Table for quantile", rq.tau[i]))
          print(coef(my.rq))
          }
         }
       }
     
     # add label
     #        print(labs)
      if(lab.pos == "bottomleft") {
        text(x = xlims[1] + diff(xlims) * .07, y = ylims[1] + diff(xlims) * .07, 
            pos = 4, cex = 1.3, col = tx.col, labels = labs)
      } else if (lab.pos == "topleft") {
        text(x = xlims[1] + diff(xlims) * .07, y = ylims[2] - diff(ylims) * .07, 
            pos = 4, cex = 1.3, col = tx.col, labels = labs)
      } else if (lab.pos == "topright") {
        text(x =  xlims[2] - diff(xlims) * .07, y = ylims[2] - diff(ylims)* .07,
            pos = 2, cex = 1.3, col = tx.col, labels = labs) 
      } else if (lab.pos == "bottomright")  {
        text(x =  xlims[2] - diff(xlims) * .07, y = ylims[1] + diff(ylims) *.07, 
            pos = 2, cex = 1.3, col = tx.col, labels = labs) 
      } 

  }

    