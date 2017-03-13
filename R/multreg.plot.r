#' A regression function created to allow for easy regression related functional performance
#' 
##  regression and multiple regression plot with one or two explanatory variables
#' @param data A input data frame, default to be NULL, 
#' @param x.add 1 to convert percent x value = 0 to 1+ x value , or 100, then *100+1, 
#' @param j threshold for logistic regression
#' @param binom if binomial function will be used 
#' @param xvar x variable, either column name or index value
#' @param x2var x variable, either column name or index value
#' @param xlog if xvar should be log transformed
#' @param xlog2 if xvar2 should be log transformed
#' @param yvar for column index of x, y values in data data frame, or column names; y could be binary
#' @param add: add multiple regression lines in single graphs espeically for logistic, regression single x variable
#' @param addcor: add principal component line between two explanatory variables and the predicted line for xvar
#' @param addref: add reference value for one variable to predict the other, e.g, add reference value for habitat(160) to predict cond
#' @param nbin: logistic regression number of bins.
#' @param LC50.plot: to add the biological response variable plot 
#' @param varlist: additional variables to include in the multiple regression plot, all these variables should be normal distributed
#' @param varpred: additional variable values for making predictions
#' @param jneg: xvar and yvar are negatively correlated, default is true
#' @param cont: if true plot contour, if not plot scatterplot 
#' @param cont2: if true plot contour, if not plot scatterplot 
#' @param model.sel  = 0 (single var) 1 (2 var), 2 (2 var interaction),3(quadratic)
#' @param rounder if x should be round to a value
#' @return Returns  
#' @keywords logistic regression, multiple regression 
#' @examples
#' set.seed(1)
#' x <- rlnorm(n) + 1:n
#' x <- (x - min(x))/max(x)
#' y <- n + 1:n + (rnorm(n))*0.2*n
#' par(mfrow=c(1,2))
#' plot(y~x)
#' scatter.plot(data = data.frame(x,y), xvar= "x", yvar = "y", add.fit = "linear", xlim = c(0,n), 
#'     ylim = range(y), log.flag = "x", x.add =100)
#' @export
#' 
multreg.plot <- function(data, yvar = "MBISQ", xvar, j = j, xvar2 = NULL, 
      jneg = TRUE, xlim = NULL, x2lim = NULL, binom = TRUE,  addcor = F, addpred = F, 
      addref = TRUE, ylim = NULL, xlab = names(data[xvar]), x2lab = names(data[xvar2]), 
      ln.col = "red", ylab= "Probability of Meeting Biological Criterion", cex.lab = 1, cex = 1,
      cex.axis = 1,  las = 1, pch = 21, col = 1, bg = 8, main = "", type = "p", cex.main = 1, 
      varlist = NULL, varpred, cont = TRUE, cont2 = FALSE,  xlog = TRUE, x2log = TRUE,  add = F, 
      rounder = 2, nbin = 41, LC50.plot = FALSE,model.sel = NULL, ...) {
 ################ Part 1 Data manipulation ################
    if(is.null(xvar2))  {           # when only one x variable
      if(is.null(data)) {          # based on independent data variable not a dataframe
         x = xvar; y = yvar
      } else {                       # data frame
         data <- data[order(data[,xvar]),]
         good <- complete.cases(data[,xvar],data[,yvar])       #   remove na values
         orig.x <- data[,xvar][good]; orig.y <- data[,yvar][good]
         x <- orig.x ; y<-orig.y        
       }
     } else {                 # if there are two x variables  or more x variables
       if(is.null(data)) {
         x = xvar; y = yvar; x2 = xvar2
       } else {    #   remove na values
         if(is.null(varlist)) {
         data <- data[order(data[,xvar]),]
         good <- complete.cases(data[,xvar], data[, xvar2], data[,yvar])
         orig.x <- data[,xvar][good]; orig.x2 <- data[,xvar2][good]; orig.y <- data[,yvar][good]
         } else {
         data <- na.omit(data[c(xvar, xvar2,yvar, varlist)])
         orig.x <- data[,xvar]; orig.x2 <- data[,xvar2]; orig.y <- data[,yvar]
#         vn <- paste("x", 3:length(varlist), sep = "")
         mycovar <- scale(data[varlist])
         covar.mean <- apply(data[varlist], 2, mean)
         covar.sd <- apply(data[varlist], 2, sd)       
         }
       x <- orig.x; x2 <- orig.x2; y <- orig.y
       }  
    }  
    print(paste("n=",length(y)))
    # log transform values
    if(xlog) {
       x <-log10(orig.x)
       }
    if(!is.null(xvar2) & x2log) {
       x2 <- log10(orig.x2)
       if(is.null(x2lim)) x2lim <- range(x2)
       else x2lim <- log10(x2lim)                       # x2lim must be at original scale
       } 

    ## Set up dummy x values for plotting
      xrange = range(x)
      if(is.null(xlim))  xlim <- xrange
      xlims <- range(xlim, xrange)

      newx = seq(xrange[1], xrange[2], length=100)
      newdata = data.frame(x = newx)
      mod <- ifelse(binom, "binomial", "gaussian")      ### determine if ordinary or logistic regression
      rtype <- ifelse(binom, "link", "response")        ## show either link or response
      if(binom)     {
       if(jneg) y <- y > j
       else  y <- y < j
       }
      crit <- ifelse(binom, 0.5, j)
      #### calculate mean for each bin
      bint <- seq(xrange[1],xrange[2], length = nbin)
      binm <- 0.5*(bint[-1] + bint[-nbin])
      binf <- cut(x, bint, include.lowest = T)
      bvals <- tapply(y, binf, mean)
   
 ############### Part 2     ####build model and predict for one x variable models
      lm0 <- glm(y ~ x, family= mod)
      lm1 <- glm(y ~ x + I(x^2), family = mod)
      if(!is.null(model.sel) & is.null(xvar2)) {
      lrm <- get(paste("lm", model.sel, sep=""))
      } else {
      if(lm0$aic > lm1$aic) {
      lrm <- lm1
      } else lrm <- lm0
      }
      predy <- predict(lrm, newdata = data.frame(x = newx), type = rtype, se.fit = T)
 #     print(head(predy))
      if(binom) {               # binomial model
        pred.y.lk <- predy$fit
        pred.up.lk <- predy$fit + qnorm(0.95) * predy$se.fit
        pred.lo.lk <- predy$fit + qnorm(0.05) * predy$se.fit
        pred.up <- exp(pred.up.lk)/(1+exp(pred.up.lk))
        pred.lo <- exp(pred.lo.lk)/(1+exp(pred.lo.lk))
        pred.y <- exp(pred.y.lk)/(1+exp(pred.y.lk))

      #### find 50% reduction from max
        max.pos <- which.max(pred.y)                     # max prob position
        mean.red <- max(pred.y)/2      # half difference
        if(max.pos == 100 | max.pos == 1) {
           mean.pos <- which.min(abs(pred.y - mean.red))  # half value position in right side data
        } else {
           mean.pos <- which.min(abs(pred.y[max.pos:100]- mean.red)) + max.pos
        }   
        LC50 <- newx[mean.pos]
     }  else {           # ordinary linear model  
        pred.y<- predy$fit
        pred.up <- predy$fit + qnorm(0.95) * predy$se.fit
        pred.lo <- predy$fit + qnorm(0.05) * predy$se.fit
      #### find j biological criterion 
        min.pos <- which.min(abs(pred.y - j))                 # find the min of the right side
        LC50 <- newx[min.pos]
        }  
      cat("One variable model", yvar, "at", j, "intercept with", xvar, "at", round(10^LC50,2), "\n")
  ############   Part 3       ############ build for two variable model   
      if(!is.null(xvar2))  {            
#### calcualte mean and sd 
      xmean <- mean(x); xsd <- sd(x);  x2mean <- mean(x2); x2sd <- sd(x2)
      x2range = range(x2)
      x2lims <- range(x2lim, x2range)
### scaled data
      x.sc <- scale(x)[,1]; x2.sc <- scale(x2)[,1]     # scale the data
      x1range.sc <- range(x.sc); x2range.sc <- range(x2.sc)
      lm0 <- glm(y~ x.sc, family = mod)   
      lm1 <- glm(y~ x.sc + x2.sc, family = mod)    
      lm2 <- glm(y~ x.sc * x2.sc, family = mod)
      lm3 <- glm(y~ x.sc + I(x.sc^2) + x2.sc + I(x2.sc^2), family = mod)
      lm4 <- glm(y~ x.sc + x2.sc + I(x2.sc^2), family = mod)
      lm5 <- glm(y~ x.sc + I(x.sc^2) + x2.sc, family = mod)
      print(c(lm0$aic, lm1$aic, lm2$aic, lm3$aic, lm4$aic, lm5$aic)) 
      min.mod <- which.min(round(c(lm0$aic, lm1$aic, lm2$aic, lm3$aic, lm4$aic, lm5$aic)))#, lm4$aic)))
      lrm <- get(paste("lm", min.mod - 1, sep=""))
      if(!is.null(model.sel)) lrm <- get(paste("lm", model.sel, sep=""))
      if(is.null(varlist))  print(summary(lrm))
      if(min.mod != "lm0") {
### build new data frame to predict 
      newx1 = seq(x1range.sc[1], x1range.sc[2], length = 100) 
      newx2 = seq(x2range.sc[1], x2range.sc[2], length = 100)      
      new.dat <- expand.grid(x.sc = newx1, x2.sc = newx2)
      pred <- predict(lrm, new.dat, type = "response",  se.fit=TRUE) 
      mean.resp <- pred$fit

 # covariance matrix between x and x2 
      xxcor <- cov(data.frame(x.sc, x2.sc))                                        
      eigenv <- eigen(xxcor)$vectors[,1]         # calculate eigenvectors, the first component of eigenvectors are extracted
      x2fit <- eigenv[2]/eigenv[1] * newx1
### alternative method is to use lmodel2 in library(lmodel2)
#      mod2 <- lmodol2(x2.sc ~x.sc )
#      x2fit <- mod2$regression.results[1,2] + mod2$regression.results[1,3]*newx1
 ### back to original scale     
      newx1.orig <- newx1 * xsd + xmean
      newx2.orig <- newx2 * x2sd + x2mean 
      newx2.fit <-  x2fit * x2sd + x2mean     
      dim(mean.resp)<- c(length(newx),length(newx2))
 ###  predict 2 variable with x2 variable fixed at mean value (0)
      new.dat2 <- data.frame(x.sc = newx1, x2.sc =rep(0,100)) 
      ynew <- predict(lrm, newdata = new.dat2, type = "response") 
       
#      xpred <- (crit - lrm$coef[1])/lrm$coef[2]
#      posit <- xpred* xsd + xmean                   good for linear regression
      posit <-  newx1.orig[which.min(abs(ynew - crit))]    # intercept in one dimension
####### predict interception between principal line and criterion line      
      pred.int <- predict(lrm, newdata = data.frame(x.sc = newx1, x2.sc = x2fit), type = "response",  se.fit=TRUE)
      min.pos <- which.min(abs(pred.int$fit - crit))  # intercept
 #     print(newx2.fit) 
      position <- c(newx1.orig[min.pos], newx2.fit[min.pos])
      pos1 <- ifelse(xlog, 10^position[1], position[1])
      pos2 <- ifelse(x2log, 10^position[2], position[2])
      cat("Principle Component Line intercept with biocriterion", yvar, "at\n",xvar, 
        xvar2, round(pos1,2), round(pos2,2),"\n")

  ########## if add reference condition for x2var
     if(is.numeric(addref)) {         
        newref <- (addref - x2mean)/x2sd
        pred.ref <- predict(lrm, newdata = data.frame(x.sc = newx1, x2.sc = rep(newref, 100)), type = "response",  se.fit=TRUE) 
        position2 <-  newx1.orig[which.min(abs(pred.ref$fit - crit))]      # 2 var intercept with reference line
        pos <- ifelse(xlog, 10^position2[1], position2[1])

  ############if add more than 2 variables
        if(!is.null(varlist)) {       # having more than 2 x variables
          scal.data <- data.frame(y, x.sc, x2.sc, mycovar)
          lrm2 <- glm(y ~ ., family = mod, data = scal.data)
          print(summary(lrm2))
          varpred.sc <- (varpred - covar.mean)/covar.sd               # scale input variable value
          addvar <- rep(varpred.sc, nrow(new.dat))
          dim(addvar) <- c(length(varpred), nrow(new.dat))
          addvar <- t(addvar)
          colnames(addvar) <- varlist
          new.dat <- data.frame(new.dat, addvar)
          position2 <- (crit - lrm2$coef[1] - sum(lrm2$coef[4:(3 + length(varlist))] * varpred.sc + lrm2$coef[3] * newref))/lrm2$coef[2]
          pos <- ifelse(xlog, 10^(position2* xsd + xmean), position2* xsd + xmean)
#        cat("Reference", xvar2, "=", addref, "intercept with biocriterion", yvar, "at\n", round(pos),"\n")

          addvar2 <- matrix(0, 100, length(varlist))                   # additional variable using mean for prediction
          colnames(addvar2) <- varlist      
          new.dat2 <- data.frame(x.sc = newx1, x2.sc =rep(0,100), addvar2) 
          ynew <- predict(lrm2, newdata = new.dat2, type = "response")  
#          posit <-  newx1.orig[which.min(abs(ynew - crit))]    # intercept in one dimension 
          posit <- ((crit - lrm2$coef[1])/lrm2$coef[2]) * xsd + xmean     
         }  
       cat("Reference", xvar2, "=", addref, "intercept with biocriterion", yvar, "at\n", (pos),"\n")   
       }
   
     LC50.2 <- posit                                       # multiplre regression intercept with j
     posit2 <- ifelse(xlog, 10^posit, posit)
     cat("Multiple regression predict", yvar, "at", j, "intercept with", xvar, "at", (posit2), "\n")
     }
    }

     ylabs <- ifelse(binom, ylab, yvar)
   #### plot model################################

     if(is.null(xvar2)) {         # if no second x variable
       if(add == FALSE) {    # if not add to another graph then plot regression line
         plot(pred.y ~ newx, type = "l", ylim = range(y, bvals, ylim, na.rm =T), 
           xlim = range(newx), axes = F, col = ln.col, main = main, xlab = xlab,
           ylab = ylabs,...)
         }  else {        # if add to another graph, then just add lines
         lines(newx, pred.y, lty = 1 , col = ln.col)
         }
         axis(2) 
         lines(newx, pred.up, lty = 2 , col = ln.col)
         lines(newx, pred.lo, lty = 2 , col = ln.col)
         points(binm, bvals, pch = pch, col = col, bg = bg, cex = cex)
         points(x, y, pch =21, cex = 0.1, col = "gray")
       }  else {           # if having more than 2 variables (x2)
        if(cont) {            # if contour then plot contour plot
           newx1.orig <<- newx1.orig;  newx2.orig <<- newx2.orig;
           ratio <<-diff(xlims)/diff(x2lims)
           mean.resp <<- mean.resp
           if(cont2) {
             image(newx1, newx2,mean.resp,
                   axes = FALSE, ann = F, xlab = xlab, ylab = x2lab)
           } else {
           eqscplot(newx1.orig, newx2.orig, ratio = diff(xlims)/diff(x2lims), type = "n",
              axes = FALSE, xlab = xlab, ylab = x2lab) 
           }
           points(x, x2, pch =21, col =1, bg = "gray", cex = cex)
          if(binom) {                # plot logistic regression
            contour(newx1.orig, newx2.orig, mean.resp, col = "black", method = "edge",
             levels = seq(0, 1, 0.1), vfont = c("sans serif", "plain"), add = TRUE,  ## original points
             axes = F, ann = FALSE, labcex=1)
           if(LC50.plot) {
            contour(newx1.orig, newx2.orig, mean.resp, col = "red", lty = 2, method = "edge",
              levels = 0.5, vfont = c("sans serif", "plain"), add = TRUE,  ## original points
              axes = F, ann = FALSE, labcex=1)
            }
          } else {                  # plot ordinary linear regression
           contour(newx1.orig, newx2.orig, mean.resp, col = "black", method = "edge",
   #          levels = seq(floor(min(mean.resp)/10)*10, ceiling(max(mean.resp)/10)*10, 5), 
             vfont = c("sans serif", "plain"), add = TRUE,  ## original points
             axes = F, ann = FALSE, labcex=1)
           if(LC50.plot) {           # add LC50 or criteria contour line
              contour(newx1.orig, newx2.orig, mean.resp, col = "red", lty = 2, method = "edge",
                 levels = j, vfont = c("sans serif", "plain"), add = TRUE,  ## original points
                 axes = F, ann = FALSE, labcex=1)
            }
          }
          if(addcor) {        # add principal component line
            lines(newx1.orig, newx2.fit, lty = 3, col = 4, lwd = 1.5)
            if(addpred) arrows(position[1], position[2], position[1], x2lims[1], lty = 3, col = 3, lwd = 3)
          }
          if(is.numeric(addref)) {          # add reference line
            abline(h = addref, lty = 2, col = 4)
            arrows(position2, addref, position2, x2lims[1], lty = 4, col = 4, lwd = 3)
          }
#        grid() 
         if(x2log) {
          max.pow2 <- ceiling(max(x2lims)); min.pow2 <- floor(min(x2lims))  ## add x range for log formation
          axis.at2 <- min.pow2:max.pow2                        # add ticks at log10 = even
          axis.lab2 <-  round(10^axis.at2, rounder)                    # add labels at log10 = even number
          axis(2, at = axis.at2, labels = axis.lab2, lwd.ticks = 1)    # major ticks with labels
          axis(2, at = log10(1:10 * rep(axis.lab2[-1]/10, each = 10)), tcl = -0.3, labels = FALSE, lwd.ticks = 0.8)  # minor ticks
          xtick2 <- axis.at2 <= max(x2lims) & axis.at2 >= x2lims[1]          # major labels
          abline(h = axis.at2, col = "lightgray", lty = "dotted")
          abline(h = log10(c(2, 5) * rep(axis.lab2[-1]/10, each = 2)),col = "lightgray", lty = "dotted")
         ### if only two or fewer major lables, then add tick labels at 2 and 5 log
          if(sum(xtick2)==2) {
           axis(2, at = log10( 5 * rep(axis.lab2[-1]/10, each = 2)), tcl = -0.4, lwd.ticks= 0.8, 
           labels = ( 5 * rep(axis.lab2[-1]/10, each = 2)))
          }
          if(sum(xtick2)<2) {
           axis(2, at = log10(c(2, 5) * rep(axis.lab2[-1]/10, each = 2)), tcl = -0.4,lwd.ticks= 0.8, 
           labels = (c(2, 5) * rep(axis.lab2[-1]/10, each = 2)))
          }
        } else {       # x2 not logged
         axis(2,las = las)
         abline(h = axTicks(2), col = "lightgray", lty = "dotted")
         }
        mtext(text = main, side = 3, line =2, cex=1)
       }  else  {   
   ##################### multiple variables and regression but only show one xvar and yvar     

         plot(newx1.orig, ynew, xlab = xlab, ylab = ylabs, ylim = range(y, ynew), axes = FALSE, type = "l", 
            col = ln.col) 
         points(x, y, pch = 20, col = "gray", cex = cex)
         axis(2)
         lines(newx, pred.y, col = "blue")                      ### original one var model
         if(binom) points(binm, bvals, pch = pch, col = col, bg = bg, cex = cex)
         if(LC50.plot) {
            abline(h = crit, lty = 2, col = 3)
            arrows(LC50.2, crit,  LC50.2, min(y), lty = 2, col = 2)     # multiple variable prediction
            arrows(LC50, crit,  LC50, min(y), col = "blue")         # one variable prediction
          }
       }
     }  
      if(xlog) {          
         max.pow <- ceiling(max(xlims)); min.pow <- floor(min(xlims))  ## add x range for log formation
         axis.at <- min.pow:max.pow                        # add ticks at log10 = even
         axis.lab <-  round(10^axis.at, rounder)                    # add labels at log10 = even number
         axis(1, at = axis.at, labels = axis.lab, lwd.ticks = 1)    # major ticks with labels
         xtick <- axis.at < max(xlims) & axis.at > xlims[1]          # major labels
         print(sum(xtick))
         abline(v = axis.at, col = "lightgray", lty = "dotted")
         abline(v = log10(c(2, 5) * rep(axis.lab[-1]/10, each = 2)),col = "lightgray", lty = "dotted")
         ### if only two or fewer major lables, then add tick labels at 2 and 5 log
         if(sum(xtick)==2) {
           axis(1, at = log10( 5 * rep(axis.lab[-1]/10, each = 2)), tcl = -0.4,lwd.ticks= 0.8, 
           labels = ( 5 * rep(axis.lab[-1]/10, each = 2)), las = las)
          }
         if(sum(xtick)<2) {
           axis(1, at = log10(c(2, 5) * rep(axis.lab[-1]/10, each = 2)), tcl = -0.4,lwd.ticks= 0.8, 
           labels = (c(2, 5) * rep(axis.lab[-1]/10, each = 2)), las = las)
          }
         if(sum(xtick) <5) {
            axis(1, at = log10(1:10 * rep(axis.lab[-1]/10, each = 10)), tcl = -0.3, labels = FALSE, lwd.ticks = 0.8)  # minor ticks
          } 
     } else    {
         axis(1,las = las)
         abline(v = axTicks(1), col = "lightgray", lty = "dotted")
     }

   box()
 }
