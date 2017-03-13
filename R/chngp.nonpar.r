#' A changing point analysis function modified from Qian, S. 's code
#' 
#' This function performed changing point analysis similar to rpart package.However, pvalues are 
#' computed in this function, as well as confidence intervals. 
#' 
#' @param data A input data frame, default to be NULL, 
#' @param bin binomial
#' @param nindex the column to define bins
#' @param num boolean variables to determine if Gausian or binomial function
#' @param xvar x variable, either column name or index value; stressor variable
#' @param yvar y variable, either column name or index value; response variable
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
#' pvalue <- chngp.nonpar(data, xvar, yvar,nindex=NULL, num = T)
#' @export

chngp.nonpar <- function(data, xvar, yvar,nindex=NULL, num = T, bin = F)
{
    if(num & bin)
        stop("data type cannot be both numeric and binary")
    data <- data [!is.na(data[,xvar])&!is.na(data[,yvar]),]
    yy <- data[,yvar]
    if (bin& max(yy)>1) stop("check data, >1 values found")
    xx <- data[,xvar]
    if (bin){
        if(is.null(data[,nindex])) nn <- rep(50, length(yy))
        else nn <- data[,nindex]
    }
    minx <- sort(unique(xx))
    m <- length(minx)                # length of unique x values
    dev <- numeric()                 # deviance storage
    if(num){
        dev [m] <- sum((yy - mean(yy))^2)    # dedevance sum(y) for all unique x values
    } else if (bin){
        nt <- sum(nn)
        p <- mean(yy)
        nt1 <- round(nt * p)
        D.1 <- ifelse(p == 0, 0, nt1*log(p))
        D.2 <- ifelse(p == 1, 0, (nt-nt1)*log(1-p))
        dev [m] <-  -2*(D.1+D.2) #deviance
    } else stop("data type must be numeric or binary")
    for(i in 1:(m-1)) {
        if(num)
            dev[i] <- sum((yy[xx <= minx[i]] - mean(yy[xx <=
                minx[i]]))^2) + sum((yy[xx > minx[i]] - mean(
                yy[xx > minx[i]]))^2)                        # dedevance sum for two groups
        else if(bin) {
            nt.left <- sum(nn[xx <= minx[i]])
            nt.right <- nt - nt.left
            p.left <- mean(yy[xx <= minx[i]])
            p.right<- mean(yy[xx >  minx[i]])
            nt.left1 <- round(nt.left * p.left)
            nt.right1<- round(nt.right* p.right)
            D.L1 <- ifelse(p.left==0, 0, nt.left1*log(p.left))
            D.L2 <- ifelse(p.left==1, 0, (nt.left- nt.left1)* log(1-p.left))
            D.R1 <- ifelse(p.right==0, 0, nt.right1* log(p.right))
            D.R2 <- ifelse(p.right==1, 0, (nt.right-nt.right1)*log(1-p.right))
            dev[i] <- -2*(D.L1 + D.L2 + D.R1 + D.R2)
        }            
    }
 ### modified so change point at the middle of two points in consistent with rpart command
    chngp <- (minx[dev == min(dev)] + minx[which(dev == min(dev)) + 1]) / 2 
    if (num){
        F.stat <- (dev[m]-min(dev))*(length(xx)-2)/(dev[m]*2)
        p.value <- round(1-pf(F.stat, 2, length(xx)-2),3)
        out <- c(chngp, p.value, mean(yy[xx <= chngp]), mean(yy[xx > chngp]))
        names(out) <- c("chngp", "p-value", "mean left","mean right")
    }else{
            nt.left <- sum(nn[xx <= chngp])
            nt.right <- nt - nt.left
            p.left <- mean(yy[xx <= chngp])
            p.right<- mean(yy[xx >  chngp])
            nt.left1 <- round(nt.left * p.left)
            nt.right1<- round(nt.right* p.right)
            D.left <- -2*(nt.left1* log(p.left) + (nt.left- nt.left1)* log(1-p.left))
            D.right<- -2*(nt.right1*log(p.right)+ (nt.right-nt.right1)*log(1-p.right))
            chi.stat <- dev[m] - D.left - D.right
            p.value <- round((1- pchisq(chi.stat, 2)),3)
        out <- c(chngp, p.value,  mean(yy[xx <= chngp]), sqrt(var(yy[xx <= chngp])),
                    mean(yy[xx >chngp]), sqrt(var(yy[xx > chngp])))    
        names(out) <- c("chngp", "p-value", "mean left", "sd left", "mean right", "sd right")
        }   
    return(out)
}
