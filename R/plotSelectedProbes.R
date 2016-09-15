# =============================================================================.
#' Plot 4C enrichment versus distance to restriction sites
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export plotSelectedProbes
#' @seealso
#'    \link{plotProbeDistanceControls},
#'    \link{enrichmentQuality}
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param a
#' @param m
#' @param dis
#' @param sel
#' @param dlim
#' @param ylab
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
plotSelectedProbes <- function(a, m, dis, sel, dlim=c(-2000, 2000), ylab="4C") {
  clr <- rgb(0,0,0,0.1)
  pch=20
  clr <- rgb((1-sel)/1.1, sel/1.1, 0, 0.1)
  plot(0, type='n', xlim=dlim, ylim=range(m, na.rm=T), xlab='probe distance', ylab=ylab)
  abline(h=0, col='lightgrey', lwd=2)
  abline(v=0, col=rgb(0.5,0,1), lwd=2)
  # points(dis, m, pch=pch, col=rgb(0,0,0,0.1), cex=0.5)
  points(dis, m, pch=pch, col=clr, cex=0.5)
  legend("topleft", "selected", fill="green", bty='n')
  legend("topright", "rejected", fill="red", bty='n')
}
