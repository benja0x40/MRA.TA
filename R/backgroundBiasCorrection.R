# =============================================================================.
#' Background bias correction
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export backgroundBiasCorrection
#' @seealso
#'    \link{backgroundBiasEstimation},
#'    \link{normalizeArrayData}
# -----------------------------------------------------------------------------.
#' @description
#' Correct for a background bias in A, M values from a tiling array.
#' The correction consists in a rotation of angle \code{-theta}, as estimated by the
#' \link{backgroundBiasEstimation} function.
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param x
#' numeric vector of A values, e.g. the average of log2 of the specific and the
#' reference signals.
#'
#' @param y
#' numeric vector of M values, e.g. the log2 ratio of the specific over the
#' reference signal.
#'
#' @param theta
#' estimated bias angle in radians.
#'
#' @param AM.scale.compensation
#' \code{logical} determining if A and M values should be rescaled prior to
#' correction by a rotation of angle \code{-theta}.
#' Scale compensation is required and activated by default.
# -----------------------------------------------------------------------------.
#' @return
#' \code{backgroundBiasCorrection} returns a \code{list} with the following
#' elements:
#' \item{x}{ corrected A value}
#' \item{y}{ corrected M value}
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
backgroundBiasCorrection <- function(x, y, theta, AM.scale.compensation = T) {
  # Scale compensation (required for A and M values)
  if(AM.scale.compensation) {
    x <- x * sqrt(2)
    y <- y * sqrt(2)/2
  }
  v <- rotate(x, y, -theta)
  v
}
