# =============================================================================.
#' Best and worst enrichment levels
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export enrichmentQuality
#' @seealso
#'    \link{plotProbeDistanceControls}
# -----------------------------------------------------------------------------.
#' @description
#' Identify probes most likely associated with best and worst enrichment levels.
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param x
#' raw Ctr signal intensity
#'
#' @param y
#' raw 4C signal intensity
#'
#' @param q.best
#'
#' @param q.worst
#'
#' @param plots
#'
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
#' x <- 2^pmin(rlnorm(100000, 1, 0.5), 10)
#' y <- 2^pmin(rlnorm(100000, 0.5, 1), 10)
#'
#' chk <- enrichmentQuality(x, y, plots=T)
#'
#' # Visualize thresholds on enrichment quality scores
#' clr <- with(
#'   chk, rgb(is.worst, is.best, 0, ifelse(is.worst | is.best, 0.5, 0.1))
#' )
#' with(chk, plot(w.score, b.score, pch='.', col=clr))
# -----------------------------------------------------------------------------.
enrichmentQuality <- function(x, y, q.best=0.015, q.worst=5E-3, plots=F) {

  n.probes <- length(x)

  A <- (log2(y) + log2(x))/2
  M <-  log2(y) - log2(x)

  qx <- rank(x)/n.probes
  qy <- rank(y)/n.probes
  qa <- rank(A)/n.probes
  qm <- rank(M)/n.probes
  # qs <- 1 - rank(abs(mean(range(A)) - A))/n.probes  # Not used

  f <- 1/2 # f = 0 => qb = qm
  qb <- (rank(qy - qx)/n.probes)^f * qm

  chk.best <- rank(qb)/n.probes > 1 - q.best

  qw <- qx * qa
  chk.worst <- rank(qw)/n.probes > 1 - q.worst
  if(plots) {
    clr <- rgb(1-qb, 1-qw, 0, 0.5)
    plot(x, y, pch='.', col=clr, main="Quality score (x,y)", log='xy')
    abline(0, 1, col="grey", lwd=2)

    clr <- rgb(chk.worst, chk.best, 0, ifelse(chk.worst | chk.best, 0.5, 0.1))
    plot(x, y, pch='.', col=clr, log='xy', main="Extrema (x,y)")
    abline(0, 1, col="grey", lwd=2)

    clr <- rgb(1-qb, 1-qw, 0, 0.5)
    plot(A, M, pch='.', col=clr, main="Quality score (A,M)")
    abline(h=0, col="grey", lwd=2)
    legend("topleft", "good", fill="green", bty='n')
    legend("topright", "poor", fill="red", bty='n')

    clr <- rgb(chk.worst, chk.best, 0, ifelse(chk.worst | chk.best, 0.5, 0.1))
    plot(A, M, pch='.', col=clr, main="Extrema (A,M)")
    abline(h=0, col="grey", lwd=2)
    legend("topleft", "best", fill="green", bty='n')
    legend("topright", "problematic", fill="red", bty='n')
  }

  list(is.best=chk.best, is.worst=chk.worst, q.best=q.best, q.worst=q.worst, b.score=qb, w.score=qw)
}
