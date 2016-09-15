# =============================================================================.
#' Normalize tiling array data
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export normalizeArrayData
#' @seealso
#'    \link{backgroundBiasEstimation},
#'    \link{backgroundBiasCorrection},
#'    \link{lowessCorrection}
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param A
#'
#' @param M
#'
#' @param smoothness
#'
#' @param epsilon
#'
#' @param nsteps
#'
#' @param name
#'
#' @param plots
#'
#' @param lowess
#'
#' @param lowess.f
#'
#' @param lowess.mad
#'
#' @param lowess.iter
#'
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
#' # Load array data and apply background bias estimation+correction, followed by
#' # lowess normalization
#'
#' data(WT.4C.Fab7.dm6)
#'
#' # 1. Combined procedure ------------------------------------------------------
#'
#' # Raw A and M values
#' A <- (log2(r1.4C$PM) + log2(r1.ct$PM))/2
#' M <- (log2(r1.4C$PM) - log2(r1.ct$PM))
#'
#' # Normalized A and M values
#' res <- normalizeArrayData(
#'   A, M, name="4C_norm", plots=TRUE
#' )
#' A <- res$A; M <- res$M
#'
#' # 2. Equivalent step by step procedure ---------------------------------------
#'
#' # Raw A and M values
#' A <- (log2(r1.4C$PM) + log2(r1.ct$PM))/2
#' M <- (log2(r1.4C$PM) - log2(r1.ct$PM))
#'
#' # Estimate background bias
#' bb.r1 <- backgroundBiasEstimation(A, M, plots = T)
#'
#' # Correct background bias
#' res <- backgroundBiasCorrection(A, M, theta=bb.r1)
#' A <- res$x; M <- res$y
#'
#' # Apply lowess normalization
#' res <- lowessCorrection(A, M, lowess.f=0.2, plots = T)
#' A <- res$x; M <- res$y
# -----------------------------------------------------------------------------.
normalizeArrayData <- function(A, M, smoothness=0.08, epsilon=0.01, nsteps=11, name="Test", plots=TRUE, lowess=T, lowess.f=0.2, lowess.mad=0, lowess.iter=5) {

  # Discard NA or infinite values
  discard <- (is.na(A) | is.na(M) | is.infinite(A) | is.infinite(M))
  ok.idx <- which(! discard)
  if(length(ok.idx)!=length(A)) {
    message("Probes associated with at least one NA in {A,M} values will not be normalized")
  }
  x <- A[ok.idx]
  y <- M[ok.idx]

  # Graphic setup
  if(plots) {
    npix <- 500
    if(lowess) {
      png(paste(name, "NormalizationReport.png",sep="_"), width=3*npix ,height=2*npix)
      layout(matrix(1:6,2,3,byrow=T))
    }
    else {
      png(paste(name, "NormalizationReport.png",sep="_"), width=2*npix ,height=2*npix)
      layout(matrix(1:4,2,2,byrow=T))
    }
    par(cex=1.4)
  }

  # Background bias estimation and correction
  theta <- backgroundBiasEstimation(
    x, y, AM.scale.compensation=T, smoothness=smoothness, epsilon=epsilon, nsteps=nsteps, plots=plots, xlab='A', ylab='M'
  )
  v <- backgroundBiasCorrection(x, y, theta, AM.scale.compensation=T)
  x <- v$x
  y <- v$y

  # MA-plot 3 (bias corrected A and M values)
  if(plots) plot(x, y ,pch='.', xlab='A', ylab='M', col=rgb(0,0,0,0.1), main="Bias corrected")

  # Lowess correction
  if(lowess) {
    v <- lowessCorrection(x, y, lowess.f=lowess.f, lowess.mad=lowess.mad, lowess.iter=lowess.iter, plots=plots)
    x <- v$x
    y <- v$y
  }

  if(plots) {
    dev.off()
  }

  A[ok.idx] <- x
  M[ok.idx] <- y

  return(list(A=A, M=M, bias=theta, ignored=which(discard)))
}
