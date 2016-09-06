# =============================================================================.
#' Compute a restriction map
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export computeRestrictionMap
#' @seealso
#'    \link{matchFragmentsToProbes},
#'    \link{matchProbesToFragments},
#'    \link{DNAStringSet}
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param genome.seq
#' \link{DNAStringSet} object containing the genome sequence
#'
#' @param enzyme.motif
#' \code{character} value defining the enzyme cut sites as a DNA sequence
#'
#' @param output.file
#' \code{character} prefix used to generate output file names.
#'
#' @param seqlist
#' indices of chromosomes to be considered (default = all).
# -----------------------------------------------------------------------------.
#' @return
#' \item{sites}{file name with the restriction sites data}
#' \item{fragments}{file name with restriction fragments data}
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
computeRestrictionMap <- function(genome.seq, enzyme.motif, output.file, seqlist=NULL) {

  out.sites     <- gsub("([.]txt)?$", "_Sites.txt", output.file, perl=T, ignore.case=T)
  out.fragments <- gsub("([.]txt)?$", "_Fragments.txt", output.file, perl=T, ignore.case=T)

  if(sum(seqlist %in% names(genome.seq))==length(seqlist) & length(seqlist)>0) {
    seqlist <- match(seqlist, names(genome.seq))
  }
  if(is.null(seqlist)) {
    seqlist <- 1:length(genome.seq)
  }

  enzyme.motif <- DNAString(enzyme.motif)

  restriction.sites     <- as.data.frame(matrix(0,0,3), stringsAsFactors=F)
  restriction.fragments <- as.data.frame(matrix(0,0,6), stringsAsFactors=F)

  for(i in seqlist) {

    seqid <- names(genome.seq)[i]
    seqln <- width(genome.seq)[i]

    motif.match <- matchPattern(enzyme.motif, genome.seq[[i]])
    n <- length(motif.match)
    message(paste(seqid, "restriction fragments\t=", n+1))

    # Restriction sites
    if(n>0) {
      restriction.sites <- rbind(
        restriction.sites, cbind(seqid, start(motif.match), end(motif.match)),
        deparse.level=0
      )
    }
    # Restriction fragments
    if(n==0) {
      rfid <- paste(seqid, ";RF", 1, sep="")
      tmp  <- cbind(rfid, seqid, 1, 1, seqln, seqln, deparse.level=0)
      restriction.fragments <- rbind(restriction.fragments, tmp, deparse.level=0)
    }
    if(n==1) {
      rfid <- paste(seqid, ";RF", 1, sep="")
      tmp <- cbind(rfid, seqid, 1, 1, start(motif.match)[1], end(motif.match)[1], deparse.level=0)
      restriction.fragments <- rbind(restriction.fragments, tmp, deparse.level=0)

      rfid <- paste(seqid, ";RF", 2, sep="")
      tmp <- cbind(rfid, seqid, start(motif.match)[n], end(motif.match)[n], seqln, seqln, deparse.level=0)
      restriction.fragments <- rbind(restriction.fragments, tmp, deparse.level=0)
    }
    if(n>1) {
      rfid <- paste(seqid, ";RF", 1, sep="")
      tmp <- cbind(rfid, seqid, 1, 1, start(motif.match)[1], end(motif.match)[1], deparse.level=0)
      restriction.fragments <- rbind(restriction.fragments, tmp, deparse.level=0)

      rfid <- paste(seqid, ";RF", 2:n, sep="")
      tmp <- cbind(rfid, seqid, start(motif.match)[1:(n-1)], end(motif.match)[1:(n-1)], start(motif.match)[2:n], end(motif.match)[2:n], deparse.level=0)
      restriction.fragments <- rbind(restriction.fragments, tmp, deparse.level=0)

      rfid <- paste(seqid, ";RF", n+1, sep="")
      tmp <- cbind(rfid, seqid, start(motif.match)[n], end(motif.match)[n], seqln, seqln, deparse.level=0)
      restriction.fragments <- rbind(restriction.fragments, tmp, deparse.level=0)
    }
  }
  colnames(restriction.sites)     <- file.formats$restriction.site$table$columns
  colnames(restriction.fragments) <- file.formats$restriction.fragment$table$columns

  write.table(format(restriction.sites, scientific=FALSE, digits=11, trim=TRUE, justify='none'), file=out.sites , sep="\t", row.names=FALSE, quote=FALSE)
  write.table(format(restriction.fragments, scientific=FALSE, digits=11, trim=TRUE, justify='none'), file=out.fragments, sep="\t", row.names=FALSE, quote=FALSE)

  n.fragments <- nrow(restriction.fragments)
  message(paste("Total computed fragments", "\t=", n.fragments))

  return(list(sites=out.sites, fragments=out.fragments))
}

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

  n.probes    <- length(probes.grg)
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
