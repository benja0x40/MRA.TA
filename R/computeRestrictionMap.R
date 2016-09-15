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
