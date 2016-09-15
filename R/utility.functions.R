# =============================================================================.
# Utility function for planar rotation (2D) centered at origin (0,0)
# -----------------------------------------------------------------------------.
# USAGE NOT DOCUMENTED
# -----------------------------------------------------------------------------.
rotate <- function(x, y, theta) {
  u <- x * cos(theta) - y * sin(theta)
  v <- x * sin(theta) + y * cos(theta)
  list(x=u, y=v)
}

# =============================================================================.
# FUNCTION line.ends
# DETAIL ----------------------------------------------------------------------.
# Utility function for visualization of the (A,M) bias slope.
# (see normalizeArrayData function below)
# -----------------------------------------------------------------------------.
# USAGE NOT DOCUMENTED
# -----------------------------------------------------------------------------.
line.ends <- function(x, y, ctr, alpha) {
  x1a <- min(x, na.rm=TRUE)
  x2a <- max(x, na.rm=TRUE)
  y1a <- ctr[2] - sin(alpha) * (ctr[1] - x1a) / cos(alpha)
  y2a <- ctr[2] + sin(alpha) * (x2a - ctr[1]) / cos(alpha)

  y1b <- min(y, na.rm=TRUE)
  y2b <- max(y, na.rm=TRUE)
  x1b <- ctr[1] - cos(alpha) * (ctr[2] - y1b) / sin(alpha)
  x2b <- ctr[1] + cos(alpha) * (y2b - ctr[2]) / sin(alpha)

  M <- c()
  if(x1b>=x1a & x1b<=x2a) M <- rbind(M, c(x1b, y1b))
  if(x2b>=x1a & x2b<=x2a) M <- rbind(M, c(x2b, y2b))
  if(y1a>=y1b & y1a<=y2b) M <- rbind(M, c(x1a, y1a))
  if(y2a>=y1b & y2a<=y2b) M <- rbind(M, c(x2a, y2a))

  if(M[1,1]<M[2,1])
    return(list(x1=M[1,1], y1=M[1,2], x2=M[2,1], y2=M[2,2]))
  else
    return(list(x1=M[2,1], y1=M[2,2], x2=M[1,1], y2=M[1,2]))
}

# =============================================================================.
#' colorstrip
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export colorstrip
#' @seealso
# -----------------------------------------------------------------------------.
#' @description
#' @details
# -----------------------------------------------------------------------------.
#' @param Pw
#'
#' @param x
#'
#' @param y
#'
#' @param width
#'
#' @param height
#'
#' @param colfun
#'
#' @param n.ticks
#'
#' @param n.bins
#'
#' @param cex
#'
#' @param \dots
#'
# -----------------------------------------------------------------------------.
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
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
