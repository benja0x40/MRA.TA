# =============================================================================.
#' Plot 4C enrichment versus distance to restriction sites
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export plotProbeDistanceControls
#' @seealso
#'    \link{plotSelectedProbes},
#'    \link{enrichmentQuality}
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param x
#'
#' @param rnk
#'
#' @param dis
#'
#' @param dlim
#'
#' @param QF
#'
#' @param n.q
#'
#' @param ylab
#'
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
plotProbeDistanceControls <- function(x, rnk, dis, dlim=c(-2000, 2000), QF=NULL, n.q=100, ylab="4C") {
  n.c <- (n.q%/%2) - 1
  a <- 0:n.c/n.c
  b <- n.c:0/n.c
  clrs <- c(
    rgb(0.5*a, 0.5*a, 1 - 0.5*a),
    rgb(0.5+0.5*a, 0.5*b^(1/10), 0.5*b)
  )
  rank.ctrl <- tapply(x, INDEX=rnk, FUN=quantile, probs=0:n.q/n.q)
  rank.ctrl <- matrix(unlist(rank.ctrl), max(rnk, na.rm=T), n.q+1, byrow=T)
  rmax <- max(which(table(rnk, useNA="no")>=3*n.q))
  rmax <- min(rmax, nrow(rank.ctrl))
  plot(0, type='n', xlim=c(1, rmax), ylim=range(x), xaxs='i', xlab='probe rank', ylab=ylab, main = "percentiles")
  grid(nx = rmax, ny=0)
  abline(h=0)
  for(i in 1:ncol(rank.ctrl)) {
    lines(1:nrow(rank.ctrl), rank.ctrl[,i], col=clrs[i])
  }
  legend("topright", c("highest","median","lowest"), fill=c("red", "grey", "blue"), bty='n')
  #boxplot(x ~ rnk, range=1, outline=T, outpch=20, outcex=0.5, outcol=rgb(0,0,0,0.1), xlab='probe rank', ylab=ylab, xlim=c(0, nrow(rank.ctrl)+1), xaxs='i')
  #abline(h=0)
  clr <- rgb(0,0,0,0.1)
  pch='.' # 46
  if(! is.null(QF)) {
    clr <- rgb(QF$is.worst, QF$is.best, 0, ifelse(QF$is.worst | QF$is.best, 0.5, 0.1))
    pch = ifelse(QF$is.worst | QF$is.best, 20, 46)
  }
  plot(0, type='n', xlim=dlim, ylim=range(x, na.rm=T), xlab='probe distance', ylab=ylab)
  abline(h=0, col='lightgrey', lwd=2)
  abline(v=0, col=rgb(0.5,0,1), lwd=2)
  points(dis, x, pch=pch, col=clr, cex=0.5)
  if(! is.null(QF)) {
    legend("topleft", "best", fill="green", bty='n')
    legend("topright", "problematic", fill="red", bty='n')
  }
}
