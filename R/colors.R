# =============================================================================.
# FUNCTION makeColors
# DETAIL -----------------------------------------------------------------------
# Compute colors from log transformed p-values with default coloring scheme
# consistent with domainogram figures in [2]. More precisely:
#   log(p value)    color gradient
#     0 to   -10	=>	white	to blue
#   -10 to   -50	=>	blue	to yellow 
#   -50 to  -500	=>	yellow	to red
#  -500 to -5000	=>	red	to black
# INPUT ------------------------------------------------------------------------
# p		vector of log transformed p-values
# thresholds	vector of numbers used to split p-values into consecutive ranges 
# colors	vector of colors defining consecutive color gradients 
# background	background color representing log(p value) = 0
# overflow	color to use for p-values exceeding the specified ranges
# OUTPUT -----------------------------------------------------------------------
# A vector of RGB colors encoded as hexadecimal character strings and reflecting
# the vector of p-values
# EXAMPLE ----------------------------------------------------------------------
# Mi <- ((cos(2*pi*(-800:800)/1600) + 1)^4)/5+3*cos(2*pi*(-800:800)/240)+2*rnorm(1601)
# Qi <- calc.Qi(Mi)
# Yi <- log(Qi)
# colors <- makeColors(Yi, thresholds=10, colors="black", overflow="red")
# -----------------------------------------------------------------------------.
makeColors <- function(p, thresholds = c(10, 50, 500, 5000), colors = c("blue","yellow","red","black"), background = "white", overflow = "white") {
  
  # Utility function, transforms a range of values into a range of colors
  p2c <- function (p, p.a, p.b, c.a, c.b) {
    i <- floor(256*(p-p.a)/(p.b-p.a))/255 
    x <- c.a[1] + i*(c.b[1]-c.a[1])
    y <- c.a[2] + i*(c.b[2]-c.a[2])
    z <- c.a[3] + i*(c.b[3]-c.a[3])
    rgb(x,y,z)
  }
  
  n.col      <- length(colors)
  colors     <- t(col2rgb(colors)/255)
  background <- as.vector(t(col2rgb(background)/255))
  
  p[p==-Inf] <- thresholds[1]
  p[is.na(p) | p==Inf] <- thresholds[length(thresholds)]
  
  p.col <- rep(rgb(1,1,1),length(p))
  
  p <- abs(p)
  thresholds <- abs(thresholds)
  
  # Colors for p-values above the least significant threshold
  x <- which(p<thresholds[1])
  if(length(x)>0) {
    p.col[x] <- p2c(p[x], 0, thresholds[1], background, as.vector(colors[1,]))
  }
  
  # Colors for intermediate p-values
  for(k in 2:n.col) {
    x <- which(p>=thresholds[k-1] & p<thresholds[k])
    if(length(x)>0) {
      p.col[x] <- p2c(p[x], thresholds[k-1], thresholds[k], as.vector(colors[k-1,]), as.vector(colors[k,]))
    }
  }
  
  # Colors for p-values below the most significant threshold
  x <- which(p>=thresholds[n.col])
  if(length(x)>0) {
    p.col[x] <- overflow
  }
  
  p.col
}

################################################################################
# FUNCTION colorstrip
# DETAIL -----------------------------------------------------------------------
# /* TODO */
# INPUT ------------------------------------------------------------------------
# /* TODO */
# OUTPUT -----------------------------------------------------------------------
# /* TODO */
# EXAMPLE ----------------------------------------------------------------------
# /* TODO */
# ------------------------------------------------------------------------------
colorstrip <- function(Pw, x, y, width, height, colfun=makeColors, n.ticks=4, n.bins=100, cex=0.5, ...) {
  
  op <- par(cex=cex)
  
  ticks <- round(min(Pw)*(0:n.ticks)/n.ticks)
  x.ticks <- x + width * (1-ticks/min(Pw))
  y.ticks <- rep(y, n.ticks)
  
  p   <- min(Pw)*(0:n.bins)/n.bins
  x.p <- x + width * (1 - (0:n.bins)/n.bins)
  y.p <- rep(y, n.bins)
  colors <- colfun(p, ...)
  
  rect(x.p[1:(n.bins-1)], y.p[1:(n.bins-1)], x.p[2:n.bins], y.p[2:n.bins]+height, col=colors, border=colors)
  
  points(c(x,x+width),c(y,y)+height,type='l')
  rect(x.ticks, y.ticks+height, x.ticks, y.ticks+height*1.3, col='black',border='black')
  text(x.ticks, y.ticks+height, ticks, pos=3)
  
  text(x, y+height/2, "log(Piw)=", pos=2)
  
  par(op)
}


