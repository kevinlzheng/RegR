#' This function calculates conditional probability values and make plots
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

cond.prob <-function(my.data, xindex, yindex, weight=NULL, goodvar = 0, goodmetric = 1, CI = TRUE,
     j,timeboot=99, xlabel=names(my.data[xindex]), ylabel = names(my.data[yindex]), cpplot = FALSE,
     log.flag="",change.pt = FALSE, rounder1=3, rounder2=1, sequencecode="",...) {

    options(scipen=10)

    n <- nrow(my.data)
    n1<- ncol(my.data)
    if(is.null(weight)) {
      weight <- n1 + 1
      my.data$wt <- rep(1,n) 
      } else my.data$wt <- my.data[, weight]

    good <- complete.cases(my.data[,xindex],my.data[,yindex], my.data[,"wt"])

    cp.data <-  my.data[good,c(xindex,yindex, "wt")]

   ################################### calculate conditional probability
    func.cond.prob <- function(data = cp.data, goodvar = goodvar, goodmetric = goodmetric, j = j) {         
      decreasing <- ifelse(goodvar == 0, TRUE, FALSE)
    levs <- sort(unique(data[,1]), decreasing = decreasing)   
    newdata = data[order(data[,1],data[,2],decreasing = decreasing),]
    n.levs <- length(levs)
    exceed.frac <- rep(NA,n.levs)
    impact.frac <- rep(NA,n.levs)          
    for(i in 1:(n.levs)){
       if(goodvar ==0) {  
         exceed.frac[i] <- sum(data[(data[,1] >= levs[i]), "wt"])
         impact.frac[i] <- ifelse(goodmetric == 1, sum(data[data[,2] <= j & data[,1] >=levs[i], "wt"]),
                    sum(data[data[, 1] >=levs[i]& data[, 2] >= j, "wt"]))
       } else {
          exceed.frac[i] <- sum(data[data[, 1] <= levs[i], "wt"])
          impact.frac[i] <- ifelse(goodmetric == 1, sum(data[data[, 2] <= j & data[,1] <=levs[i], "wt"]),
                   sum(data[data[, 1] <=levs[i] & data[, 2] >= j , "wt"]))
       }
      }
     cond.prob <- impact.frac/exceed.frac
     paired<-data.frame(levs,cond.prob)
    return(paired)
    }
   ##################################  conditional probability calculation ends 
    
    cp <- func.cond.prob(data = cp.data, goodvar, goodmetric, j)
    print(head(cp))
      w <- ifelse(goodmetric == 1,"<=", ">=")
      if(goodvar == 0) xx <- range(cp.data[, 1])
      else xx = rev(range(cp.data[, 1]))
 #   if(CI == FALSE)   {
      plot(cp$levs, cp$cond.prob, axes = F, type = "p", xlab = xlabel, col="darkgreen",
               ylab = paste("Prob. of",ylabel,w, round(j,rounder2), sep=" "), xlim = xx,
               log = log.flag, xpd=F)
       axis(1, at=axTicks(1), labels=round(axTicks(1),rounder1) )
       axis(2, at=axTicks(2), labels=round(axTicks(2),2))
       box()
#    } else   {
    if(CI == TRUE) {                 
    for(index in 1:timeboot){
      isamp <- sample(1:nrow(cp.data), replace = TRUE, prob = cp.data[, "wt"])
      data.resamp <- cp.data[isamp, ]
      cond.prob <- func.cond.prob(data = data.resamp, goodvar, goodmetric, j)
      if(index == 1) cp.boot <- cp
      cp.boot <-merge(cp.boot,cond.prob, by.x="levs", by.y="levs",all=T)
      }
      print(head(cp.boot))
    lowCI<-apply(cp.boot[,-1],1,quantile,0.05, na.rm = TRUE)
    highCI<-apply(cp.boot[,-1],1,quantile,0.95, na.rm = TRUE)
    middle<-apply(cp.boot[,-1],1,quantile,0.5, na.rm = TRUE)
    allCI<-data.frame(cbind(cp.boot$levs, lowCI,middle,highCI))
     
    lines(cp.boot$levs, highCI,lty="dashed",lwd=1.5,col="gray50")
    lines(cp.boot$levs, lowCI, lty="dashed",lwd=1.5,col="gray50")
    text(x=min(cp.boot$levs)*1.1, y = max(middle)*0.95, cex=1.5, col=2, labels=sequencecode) 

 #   if(change.pt)   {
#    source("C:/Users/lei.zheng/Desktop/R_Function/func.change.ci.r")
#    ci <- func.change.ci(cp, 1, 2)                          # calculate the changepoint
 #   abline(v=ci[1],lty="dashed",col="gray50", lwd=1.5)            # add the change point line
#    abline(v=ci[2],lty="solid",lwd=1.5)                           # add the change point line
 #   abline(v=ci[3],lty="dashed",col="gray50", lwd=1.5)            # add the change point line
#    mtext(paste(round(ci[1],rounder1), round(ci[2],rounder1),round(ci[3],rounder1),sep="~"), side=3, line=0.3)
#    }
   }    
  }
