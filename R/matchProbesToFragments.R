# =============================================================================.
#' Match micro array probes to restriction fragments
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export matchProbesToFragments
#' @seealso
#'    \link{matchFragmentsToProbes},
#'    \link{computeRestrictionMap}
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param probes.grg
#'
#' @param fragments.grg
#'
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
matchProbesToFragments <- function(probes.grg, fragments.grg) {
  # ---------------------------------------------------------------------------.
  match.fragments <- function(ovl, probes.middle, regions.start, regions.end, regions.middle) {

    res <- as.data.frame(matrix(NA, queryLength(ovl), 6))
    names(res) <- c("rfi", "f.nbr", "site", "side", "rnk", "dst")

    qry <- queryHits(ovl)
    sbj <- subjectHits(ovl)
    chk <- probes.middle[qry]>regions.start[sbj] & probes.middle[qry]<=regions.end[sbj]
    ovl <- ovl[chk]; qry <- qry[chk]; sbj <- sbj[chk];

    res$rfi[qry]  <- sbj
    res$f.nbr <- countQueryHits(ovl)

    side <- ifelse(probes.middle[qry]>regions.middle[sbj],2,1)
    idx  <- 1:length(sbj) + length(sbj)*(side-1)
    res$site[qry] <- cbind(regions.start[sbj], regions.end[sbj])[idx]
    res$side[qry] <- ifelse(side==1,5,3)

    chk <- ! duplicated(paste(sbj, side))
    idx <- cumsum(chk)
    nbr <- diff(c(which(chk), length(idx)+1))[idx]

    dst <- probes.middle[qry] - res$site[qry]
    rnk <- unlist(by(dst, INDICES=idx, FUN=order))
    chk <- dst<0
    rnk[chk] <- nbr[chk] - rnk [chk] + 1

    res$dst[qry] <- dst
    res$rnk[qry] <- rnk

    res
  }
  # ---------------------------------------------------------------------------.
  message("Matching probes to restriction fragments")

  probes.grg    <- sort(probes.grg)
  fragments.grg <- sort(fragments.grg)

  probes.grg$MIDDLE    <- (start(probes.grg)+end(probes.grg))/2
  fragments.grg$MIDDLE <- (start(fragments.grg)+end(fragments.grg))/2

  ovl <- findOverlaps(probes.grg, fragments.grg, maxgap=0, minoverlap=1, type="any")

  res <- match.fragments(ovl, probes.grg$MIDDLE, start(fragments.grg), end(fragments.grg), fragments.grg$MIDDLE)

  probes.grg$RF_ID    <- names(fragments.grg)[res$rfi] # fragment id
  probes.grg$RF_START <- start(fragments.grg)[res$rfi] # fragment id
  probes.grg$RF_MID   <- end(fragments.grg)[res$rfi] # fragment id
  probes.grg$RF_END   <- end(fragments.grg)[res$rfi] # fragment id
  probes.grg$RF_LEN   <- width(fragments.grg)[res$rfi] # fragment length
  probes.grg$RF_SITE  <- res$site # nearest restriction site
  probes.grg$RF_SIDE  <- res$side # probe is within the 5' or 3' end of the fragment
  probes.grg$RF_RANK  <- res$rnk  # rank of probe positions according to nearest restriction site
  probes.grg$RF_DIST  <- res$dst  # distance to nearest restriction site

  probes.grg
}
