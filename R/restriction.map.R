
# =============================================================================.
#' Count the number of accepted probes per restriction fragment
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export countProbesPerFragment
#' @seealso
#'    \link{matchProbesToFragments},
#'    \link{matchFragmentsToProbes},
#'    \link{computeRestrictionMap}
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param probes.grg
#'
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
countProbesPerFragment <- function(probes.grg) {
  n.probes <- length(probes.grg)
  apnbr <- by(1:n.probes, INDICES = probes.grg$RF_ID, FUN = length)
  idx <- match(probes.grg$RF_ID, names(apnbr))
  probes.grg$RF_APNBR <- apnbr[idx]
  probes.grg
}

# =============================================================================.
#' Combine probe measurements by restriction fragments
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export combineByFragment
#' @seealso
#'    \link{}
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param x
#'
#' @param probes.grg
#'
#' @param FUN
#'
#' @param \dots
#'
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
combineByFragment <- function(x, probes.grg, FUN = mean, ...) {

  fs <- paste(probes.grg$RF_ID, probes.grg$RF_SIDE, sep="_")
  s <- unlist(by(probes.grg$RF_SITE, INDICES=fs, FUN=function(k) { unique(k) }))
  y <- unlist(by(x, INDICES=fs, FUN=FUN, ...))
  data.frame(
    RF_ID = gsub("^(.+)_([53])$", "\\1", names(y), perl=T),
    SIDE  = gsub("^(.+)_([53])$", "\\2", names(y), perl=T),
    SITE  = as.vector(s),
    VALUE = as.vector(y),
    stringsAsFactors=F
  )
}
