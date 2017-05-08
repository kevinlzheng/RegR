#' A three dimensional plot function
#' 
#' The function persp is the base graphics function for creating wireframe surface plots.
#' adopted from http://www.r-bloggers.com/creating-surface-plots/
#' The persp function requires a list of x and y values covering the grid of vertical values
#'  which is specified as the z variable. The heights for the display are specified as a table
#'  of values which we saved previously as the object z during the calculations when the local
#'  trend surface model was fitted to the data.
#' 
#' @param x  
#' @param y a boolean variable deciding if the box should be horizontal 
#' @param z e if an outline 
#' @param xlab
#' @param ylab
#' @param zlab
#' @param phi 
#' @param theta The function arguments phi and theta are used to rotate the viewing angle of the surface
#' @param border border of lines
#' @param box should hte bounding box for the surface be displayed
#' @example 
#' set.seed(200)
#' n = 20
#' x <-  1:n  # rlnorm(n) +
#' x <- (x - min(x))/max(x)
#' x <- rep(x, n)
#' y <- n + 1:n # + (rnorm(n))*0.2*n
#' y <- rep(y, each = n)
#' z <- rowSums(expand.grid(x = rnorm(n)+1:n*2, y = 1:n + (rnorm(n, 1, 3))) )  # c("gray1", "gray99")
#' scheme = c("purple4", "blue4","blue1","green3", "yellowgreen","yellow2","sandybrown","red","darkred")
#' surf3d(x,y,z, type ="cm")
surf3d <- function(x, y, z, xlab = "X Coordinate (feet)", ylab = "Y Coordinate (yyy)", border = NA,
                   zlab ="Elevation", main ="Surface elevation data", phi=45, theta =45,ncut = 100,
                   scheme = c("purple4", "blue4","blue1","green3", "yellowgreen","yellow2","sandybrown","red","darkred"),
                   type ="rgb") {
 #   ncut <- 11
    newx <- seq(min(x), max(x), length = ncut)
    newy <- seq(min(y), max(y), length = ncut)
    lm1 <- loess(z ~ x * y, degree = 2, span = .25)
    newz <- predict(lm1, newdata = expand.grid(x = newx, y = newy))
    rgb.palette <- colorRampPalette(scheme, space = "rgb")

    # Compute the z-value at the facet centres
    zfacet <- newz[-1, -1] + newz[-1, -ncut] + newz[-ncut, -1] + newz[-ncut, -ncut]
    # Recode facet z-values into color indices
    # facetcol <- cut(zfacet, ncut)
    brks <- seq(min(zfacet), max(zfacet), length = ncut)    
    ranks <- matrix(as.numeric(cut(zfacet, brks, include.lowest =T)),
                   nrow = ncut-1, ncol = ncut-1, byrow = FALSE)
    rgb.palette <- colorRampPalette(scheme, space = "rgb")  
    col.val <- as.numeric(cut(zfacet, brks, include.lowest =T))
    val.col <- switch(type, 
                      rgb = rgb.palette(ncut),
                      topo = topo.colors(ncut),
                      heat = heat.colors(ncut),
                      rain = rainbow(ncut),
                      cm = cm.colors(ncut),  
                      ter = terrain.colors(ncut)) 

    cols <- matrix(val.col[col.val], nrow = ncut-1, ncol = ncut -1, byrow = FALSE)
    persp(newx, newy, (newz), phi = phi, theta = -30, border = border, 
        xlab = xlab, ylab = ylab, zlab =zlab, main = main,
        col = cols, ...)
    }
