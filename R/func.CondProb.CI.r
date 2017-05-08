#' This function calculates confidence interval of conditional probability values and plots
#' 
#' It also uses bootstrapping techniques to calculate confidence interval of change points 
#' 
#' @param my.data data frame 
#' @param xindex xvariable at column name or index 
#' @param yindex xvariable at column name or index 
#' @param goodmetric is a flag for good or bad metrics, e.g., EPT metric is good(1), HBI is bad (0).
#' @param weight is also optional, when there is a site weight column in mydata, enter the column number
#' @param cpplot is a optional flag to determine if conditional prob should be plotted conditional prob with confidence intervals should be plotted.
#' @param change.pt is an option flag to choose if changepoint analysis should be performed
#' @param goodvar is a flag for good or bad variables, goodvar is 1 if a predictor increase a good metric, e.g., habitat score
#' @param rounder1 is number of decimal for change points displaying in graphs
#' @param rounder2 is biocriteria rounder
#' @param j x is stressor, y is biological response, j is the criteria
#' @examples
#' 
#' @export

func.CondProb.CI<-function(my.data, xindex, yindex, weight=NULL, goodvar = FALSE, goodmetric = TRUE, CI = TRUE,
     j,timeboot=99, xlabel=names(my.data[xindex]), ylabel = names(my.data[yindex]), cpplot = FALSE,
     log.flag="",change.pt = FALSE, rounder1=3, rounder2=1, sequencecode="",...) {

    options(scipen=10)

    n <- nrow(my.data)
    n1<- ncol(my.data)
    if(is.null(weight)) {
      weight <- rep(1,n)
      my.data$weight = weight 
      }
    good <- complete.cases(my.data[,xindex],my.data[,yindex], y.data[,weight])
    cp.data <-  my.data[good,c(xindex,yindex,weight)]
   ################################### calculate conditional probability
    func.cond.prob <- function(data = cp.data, goodvar = goodvar, goodmetric = goodmetric, j = j) {         
      decreasing <- ifelse(goodvar == 0, TRUE, FALSE)
    levs <- sort(unique(my.data[,1]), decreasing = decreasing)   
    newdata = cp.data[order(cp.data[,1],cp.data[,2],decreasing = decreasing),]
    n.levs <- length(levs)
    exceed.frac <- rep(NA,n.levs)
    impact.frac <- rep(NA,n.levs)          
    for(i in 1:(n.levs)){
       if(goodvar ==0) {  
         exceed.frac[i] <- sum(cp.data[(cp.data[,xindex] >= levs[i]), weight])
         impact.frac[i] <- ifelse(goodmetric == 1, sum(cp.data[cp.data[,yindex] <= j & cp.data[,xindex] >=levs[i],weight]),
                    sum(cp.data[cp.data[,xindex] >=levs[i]& cp.data[,yindex] >= j, weight]))
       } else {
          exceed.frac[i] <- sum(cp.data[cp.data[,xindex] <= levs[i],weight])
          impact.frac[i] <- ifelse(goodmetric == 1, sum(my.data[my.data[,yindex] <= j & my.data[,xindex] <=levs[i],weight]),
                   sum(my.data[my.data[,xindex] <=levs[i] & my.data[,yindex] >= j , weight]))
       }
      }
     cond.prob <- impact.frac/exceed.frac
     paired<-data.frame(levs,cond.prob)
    return(paired)
    }
   ##################################  conditional probability calculation ends 
    cp <- func.cond.prob(data = cp.data, goodvar, goodmetric, j)
   
      ifelse(goodmetric == 1,w = "<=", w = ">=")
      
      if(goodvar == 0) xx <- range(cp.data[, 1])
      else xx = rev(range(cp.data[, 1]))
      plot(levs, cond.prob, axes = F, type = "p", xlab = xlabel, col="darkgreen",
               ylab = paste("Prob. of",ylabel,w, round(j,rounder2), sep=" "), xlim = xx, log = log.flag, xpd=F)
       axis(1, at=axTicks(1), labels=round(axTicks(1),rounder1) )
       axis(2,at=axTicks(2), labels=round(axTicks(2),2))
      box()
    
      if(cpplot==TRUE)   {
        func.cond.prob(my.data, xindex, yindex, weight, goodvar, goodmetric, j,
                  xlabel, ylabel, cpplot=cpplot, rounder1=rounder1, rounder2=rounder2, 
                       log=log.flag,sequencecode="")
        } else {
        cond.prob.boot<-func.cond.prob(my.data, xindex, yindex, weight, goodvar, goodmetric, j,
                  xlabel, ylabel, cpplot=cpplot, log=log.flag)
      }           
    for(index in 1:timeboot){
      isamp <- sample(1:nrow(my.data), replace = TRUE, prob=(!is.na(my.data[, weight])))
      my.data.resamp <- my.data[isamp, ]
      cond.prob <- func.cond.prob(my.data.resamp, xindex, yindex, weight, goodvar = goodvar, 
          goodmetric = goodmetric, j, cpplot = FALSE,log = log.flag)
      cond.prob.boot<-merge(cond.prob.boot,cond.prob, by.x="levs", by.y="levs",all=T)
      }

    lowCI<-apply(cond.prob.boot[,-1],1,quantile,0.05, na.rm = TRUE)
    highCI<-apply(cond.prob.boot[,-1],1,quantile,0.95, na.rm = TRUE)
    middle<-apply(cond.prob.boot[,-1],1,quantile,0.5, na.rm = TRUE)
    allCI<-data.frame(cbind(cond.prob.boot$levs, lowCI,middle,highCI))
#    myref<-deparse(substitute(my.data))
    ifelse(goodmetric == 1,w<-"<=", w<-">=")
    xx<-range(my.data[,xindex], na.rm=TRUE)
    ifelse(goodvar == 0,xx, xx<-c(xx[2],xx[1]))

    plot(cond.prob.boot$levs, middle, axes = F, type = "p", xlab ="", xlim = xx, 
         col="darkgreen", ylab = "",log=log.flag, xpd=F)   
       axis(1, at=axTicks(1), labels=round(axTicks(1),rounder1) )
       axis(2, at=axTicks(2), labels=round(axTicks(2),2))
     box() 
       
    lines(cond.prob.boot$levs, highCI,lty="dashed",lwd=1.5,col="gray50")
    lines(cond.prob.boot$levs, lowCI, lty="dashed",lwd=1.5,col="gray50")
    mtext(paste("Prob. of",ylabel,w, round(j, rounder2)), side = 2, line = 2.3)
    text(x=min(cond.prob.boot$levs)*1.1, y=max(middle)*0.95, cex=1.5, col=2, labels=sequencecode) 
    mtext(xlabel, side = 1, line = 2.3)

#    if(change.pt)   {
#    ci<-func.change.ci(allCI,1,3)                          # calculate the changepoint
#    abline(v=ci[1],lty="dashed",col="gray50", lwd=1.5)            # add the change point line
#    abline(v=ci[2],lty="solid",lwd=1.5)                           # add the change point line
#    abline(v=ci[3],lty="dashed",col="gray50", lwd=1.5)            # add the change point line
#    mtext(paste(round(ci[1],rounder1), round(ci[2],rounder1),round(ci[3],rounder1),sep="~"), side=3, line=0.3)
#    }
  
  }
