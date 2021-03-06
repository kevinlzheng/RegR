% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/box.plot.r
\name{box.plot}
\alias{box.plot}
\title{A box plot function showing pretty box plot and jitter data points}
\usage{
box.plot(x, horizontal = TRUE, outline = FALSE, lty = 1, staplewex = 0,
  line = 6, boxwex = 0.8, boxlwd = 1, medlwd = 1, xaxt = "n",
  yaxt = "n", notch = FALSE, mar = c(5, 8, 3, 1),
  xlab = "Assigned Probability \%", ylab = "Phrase",
  main = "Perceptions of Probability", cex.main = 1, names = NULL)
}
\arguments{
\item{x}{A data frame or matrix showing with cross tabled format,}

\item{horizontal}{a boolean variable deciding if the box should be horizontal}

\item{outline}{if an outline}

\item{line}{when mar changes, ylab should be changed to by line = 3}

\item{mar}{margin for plot area}

\item{main}{for title}

\item{names}{group names, can be entered manually}
}
\description{
This function box plot similar to ggplot2 function
modified from https://github.com/zonination/perceptions/blob/master/percept.R)
}
\examples{
set.seed(100)
n <- rnorm(11, 2, 11)
x <-  (matrix(rnorm(15*11, 3, 7), nrow=15, ncol=11) *2)
x2 <- sapply(1:11, function(i) x[,i]+n[i]) 
names <- paste("Var", 1:11, sep =".")
box.plot (x2, horizontal = TRUE, outline = FALSE, lty=1, staplewex=0, 
         boxwex=0.8, boxlwd=1, medlwd=1,xaxt="n", yaxt="n", notch = FALSE,
         mar=c(5,5,3,1), xlab ="Assigned Probability \%", ylab = "Phrase",
         main = "Perceptions of Probability", names = names) 
}
