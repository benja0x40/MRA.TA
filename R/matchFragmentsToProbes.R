# =============================================================================.
#' Match restriction fragments to microarray probes
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export matchFragmentsToProbes
#' @seealso
#'    \link{matchProbesToFragments},
#'    \link{computeRestrictionMap}
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param fragments.grg
#'
#' @param probes.grg
#'
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
matchFragmentsToProbes <- function(fragments.grg, probes.grg) {
  # ---------------------------------------------------------------------------.
  match.probes <- function(ovl, probes.middle, regions.start, regions.end) {

    res <- as.data.frame(matrix(NA, subjectLength(ovl), 3))
    names(res) <- c("nbr", "first", "last")

    qry <- queryHits(ovl)
    sbj <- subjectHits(ovl)
    chk <- probes.middle[qry]>regions.start[sbj] & probes.middle[qry]<=regions.end[sbj]

    res$nbr <- countSubjectHits(ovl[chk])
    idx <- by(qry[chk], INDICES=sbj[chk], FUN=min)
    res$first[as.numeric(names(idx))] <- as.vector(idx)
    idx <- by(qry[chk], INDICES=sbj[chk], FUN=max)
    res$last[as.numeric(names(idx))] <- as.vector(idx)

    res
  }
  # ---------------------------------------------------------------------------.
  probes.grg    <- sort(probes.grg)
  fragments.grg <- sort(fragments.grg)

  probes.grg$MIDDLE    <- (start(probes.grg)+end(probes.grg))/2
  fragments.grg$MIDDLE <- (start(fragments.grg)+end(fragments.grg))/2

  fragments.grg$PROBES_NBR <- 0  # Number of matching probes
  fragments.grg$START_IDX  <- NA # Matching probe that is the closest to 5' restriction site
  fragments.grg$RS5_IDX    <- NA # 5' matching probe that is the closest to fragment center
  fragments.grg$RS3_IDX    <- NA # 3' matching probe that is the closest to fragment center
  fragments.grg$END_IDX    <- NA # Matching probe that is the closest to 3' restriction site
  fragments.grg$RS5_NBR    <- 0  # Number of matching probes from 5' site to center
  fragments.grg$RS3_NBR    <- 0  # Number of matching probes from center to 3' site

  message("Matching restriction fragments to probes")

  ovl <- findOverlaps(probes.grg, fragments.grg, maxgap=0, minoverlap=1, type="any")

  res <- match.probes(ovl, probes.grg$MIDDLE, start(fragments.grg), end(fragments.grg))
  fragments.grg$PROBES_NBR <- res$nbr
  fragments.grg$START_IDX  <- res$first
  fragments.grg$END_IDX    <- res$last

  res <- match.probes(ovl, probes.grg$MIDDLE, start(fragments.grg), fragments.grg$MIDDLE)
  fragments.grg$RS5_NBR <- res$nbr
  fragments.grg$RS5_IDX <- res$last

  res <- match.probes(ovl, probes.grg$MIDDLE, fragments.grg$MIDDLE, end(fragments.grg))
  fragments.grg$RS3_NBR <- res$nbr
  fragments.grg$RS3_IDX <- res$first

  fragments.grg
}
