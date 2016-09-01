# =============================================================================.
#'
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export
#' @seealso
#'    \link{},
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
localopt.c <- function(Yi, name, wmax = 0, wmin = 0, gamma = 1, Tw = 0) {

  # Count number of measurements
  n.probes <- length(Yi)

  # Default minimum domain size = 3
  wmin <- max(2, wmin-1)

  # Default maximum window size = n.probes
  if(wmax < wmin + 1)
    wmax <- n.probes
  wmax <- min(wmax, n.probes)

  # Default thresholds = 0
  if(length(Tw)==1)
    Tw <- rep(Tw, n.probes)

  message("[ localopt ] name = ", name, "\tn.probes = ", n.probes, "\twmax = ", wmax, "\twmin = ", wmin+1, "\tgamma = ", gamma)

  # ---------------------------------------------------------------------------.
  # Output file
  output.name <- paste(name, "_RawSegments.txt", sep="")

  # ---------------------------------------------------------------------------.
  # Optimization
  Raw.segments <- .C(
    "localopt_fisher",
    as.double(Yi),        # INPUT: sequence of measurements
    as.integer(n.probes), # INPUT: number of measurements
    as.character(name),   # INPUT: name used for the output file
    as.integer(wmax),     # INPUT: maximal scale to explore as number of measurements
    as.integer(wmin),     # INPUT: minimal domain size as number of measurements
    as.double(gamma),	    # INPUT: p-value improvement factor
    as.double(Tw),        # INPUT: p-value threshold for each scale w
    Pw.min = double(n.probes),		 # OUTPUT: best p-value for each scale w
    n.segments = integer(1)  # OUTPUT: number of optimal domains
  )

  return(list(Pw.min=Raw.segments$Pw.min, n.segments=Raw.segments$n.segments, file=output.name))
}
