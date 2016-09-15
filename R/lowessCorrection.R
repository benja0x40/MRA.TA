# =============================================================================.
# *** WARNING/TODO *** #  Cleanup the returned A value (and list names?)
# =============================================================================.
#' LOWESS based normalization
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export lowessCorrection
#' @seealso
#'    \link{normalizeArrayData}
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param x
#'
#' @param y
#'
#' @param lowess.f
#'
#' @param lowess.mad
#'
#' @param lowess.iter
#'
#' @param plots
#'
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
lowessCorrection <- function(x, y, lowess.f=0.3, lowess.mad=0, lowess.iter=5, plots=F) {
  A <- x
  M <- y

  if(lowess.mad>0) {
    j <- which(abs(y-median(y))<lowess.mad*mad(y))
    normalization <- lowess(x[j], y[j], f=lowess.f, iter=lowess.iter)
    o <- order(A)
    yn <- approx(normalization$x, normalization$y, xout=A[o], rule=2, ties='ordered')$y
    y <- M - yn[order(o)]
    # x <- A + yn[order(o)]/2
  }
  else {
    j <- 1:length(x)
    normalization <- lowess(x, y, f=lowess.f, iter=lowess.iter)
    o <- order(A)
    yn <- normalization$y
    y <- M - yn[order(o)]
    # x <- A + yn[order(o)]/2
  }
  if(plots) {
    title <- paste("Lowess (f =", lowess.f)
    if(lowess.mad==0) title <- paste(title,")",sep="")
    else title <- paste(title, ", ", lowess.mad," x mad)", sep="")
    # MA-plot 4 (overlayed with a visualization of the lowess)
    plot(A, M ,pch='.', col=rgb(0,0,0,0.1), main=title)
    if(lowess.mad>0)
      points(A[j], M[j] ,pch='.', col=rgb(0,0.25,0.75,0.05))
    lines(A[o], yn, col=rgb(0.5,0,1,1))
    # MA-plot 5 (overlayed with a visualization of the lowess)
    plot(x, y ,pch='.', col=rgb(0,0,0,0.1), main="Lowess corrected", xlab='A', ylab='M')
    lines(c(min(x),max(x)),c(0,0),col=rgb(0.3,0.3,0.3,0.7),lwd=1.5)
  }
  list(x=x, y=y)
}
