# =============================================================================.
#' Multi-resolution segmentation
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export segmentation
#' @seealso
#'    \link{calc.Qi},
#'    \link{enrichmentScore},
#'    \link{domainogram},
#'    \link{plotOptimalSegments},
#'    \link{plotDomains}
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param Yi
#' vector of statistical scores produced by the \link{enrichmentScore} function.
#'
#' @param name
#' prefix used to generate output file names.
#'
#' @param wmax
#' maximal window size to be scanned in number of Yi values.
#' The default is \code{length(Yi)}.
#'
#' @param wmin
#' minimum size of segmented domains in number of Yi values.
#' The default is 3 consecutive Yi values.
#'
#' @param gamma
#' stringency factor ranging from ]0;1].
#' Default value is 1.0 and lower values correspond to increased stringency.
#'
#' @param Tw
#' log value ranging from ]\code{-Inf};0] and indicating the minimal statistical
#' scores (e.g. pseudo p-values from combined rank-based scores) accepted for
#' segmented domains.
#' Default value is 0 and equivalent to an absence of cutoff while lower values
#' correspond to more stringent cutoffs.
# -----------------------------------------------------------------------------.
#' @return
#' \code{segmentation} returns a list object with the following attributes:
#' \item{time.optimization}{computation time for the optimization procedure}
#' \item{time.MRT.Analysis}{computation time for the domain fusion procedure}
#' \item{n.segments}{total number of locally optimal segments after optimization}
#' \item{n.domains}{total number of domains after domain fusion}
#' \item{n.maxresolution}{number of maximum resolution domains}
#' \item{n.maxscale}{number of maximum scale domains}
#' \item{file.segments}{result file for the optimization procedure}
#' \item{file.domains}{result file for the domain fusion procedure}
#' \item{file.maxresolution}{result file for maximum resolution domains}
#' \item{file.maxscale}{result file for maximum scale domains}
#' \item{Pw.min}{best statistical score (pseudo p-values from combined rank-based scores) for each scanned window size}
# -----------------------------------------------------------------------------.
#' @section
#' Output files: \itemize{
#' \item Common structure
#'
#' All output files are tab delimited text files including two header lines.
#' \preformatted{ line 1 = analysis parameters line 2 = column names following
#' lines = data table }
#'
#' \item \code{file.segments}
#'
#' This file contains the optimal segment data, resulting from the the
#' optimization procedure, and defined by columns (i, w, Piw) as follows:
#' \preformatted{ i = index of the last Yi value within the optimal segment w =
#' size in number of consecutive Yi values Piw = statistical score }
#'
#' \item \code{file.domains}, \code{file.maxresolution}, \code{file.maxscale}
#'
#' These files contain the segmented domain data, resulting from the domain
#' fusion procedure, and defined by columns (id, container, start, end, wmin,
#' wmax, P, i, w) as follows: \preformatted{ id = unique identifier for each
#' domain container = identifier of parent domain, 0 meaning no parent start =
#' index of the first Yi value within domain end = index of the last Yi value
#' within domain wmin = minimum size of included optimal segments, in number of
#' Yi values wmax = maximum size of included optimal segments, in number of Yi
#' values P = statistical score of the locally optimal segment i = index of the
#' last Yi value within the locally optimal segment w = size of the locally
#' optimal segment in number of consecutive Yi values }
#' }
# -----------------------------------------------------------------------------.
#' @examples
#' # Simulate enrichment signal
#' n <- 2000
#' Mi <- rep(0, n)
#' Mi <- Mi + dnorm(1:n, 2.5*n/20,  n/40) + dnorm(1:n, 4*n/20,  50)
#' Mi <- Mi + 4 * dnorm(1:n, 5*n/10, n/10)
#' Mi <- Mi + dnorm(1:n, 16*n/20,  n/40) + dnorm(1:n, 17.5*n/20,  50)
#' Mi <- (Mi/max(Mi))^4 + rnorm(n)/4
#'
#' # Compute enrichment scores
#' Yi <- enrichmentScore(Mi)
#'
#' # Multi-resolution segmentation
#' seg.c <- segmentation(Yi, name="MRA_demo", wmin=20)
#'
#' # Load segmentation results
#' opts <- read.delim(seg.c$file.segments, stringsAsFactors=F, skip=1)
#' doms <- read.delim(seg.c$file.domains, stringsAsFactors=F)
#' doms.mr <- read.delim(seg.c$file.maxresolution, stringsAsFactors=F)
#' doms.ms <- read.delim(seg.c$file.maxscale, stringsAsFactors=F)
#'
#' # Visualization coordinates
#' x.s <- 1:n - 0.5
#' x.e <- 1:n + 0.5
#' w2y <- wSize2yAxis(n, logscale=T)
#'
#' layout(matrix(1:2, 2, 1), heights=c(3,1)/4)
#' par(mar=c(3, 4, 1, 2)) # default bottom, left, top, right = c(5, 4, 4, 2)
#'
#' # Plot domainogram
#' domainogram(Yi, x.s, x.e, w2y)
#' plot(Mi, type='l', xaxs = 'i')
#'
#' # Visualize segmentation results
#' plotOptimalSegments(opts, x.s, x.e, w2y, col="black")
#' plot(Mi, type='l', xaxs = 'i')
#'
#' # Visualize multi-resolution domains
#' plotDomains(doms, x.s, x.e, w2y, col=rgb(0,0,0,0.5), border=rgb(0,0,0,0))
#' # Visualize max. resolution and max. scale domains
#' plotDomains(doms.mr, x.s, x.e, w2y, col=rgb(0,1,0,0.5), border=rgb(0,0,0,0), lwd=2, lty=1, add=T)
#' plotDomains(doms.ms, x.s, x.e, w2y, col=rgb(1,0,0,0.5), border=rgb(0,0,0,0), lwd=2, lty=1, add=T)
#' legend("topright", c("Max. resolution", "Max. scale", "Both"), fill=c("green", "red", "chocolate"), bty='n')
#' plot(Mi, type='l', xaxs = 'i')
# -----------------------------------------------------------------------------.
segmentation <- function(Yi, name, wmax=0, wmin=3, gamma=1, Tw=0) {

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

  # Check basic parameter consistency
  if(n.probes <= wmin)
    stop(paste("number of measurements Yi should be greater or equal to", wmin + 1))
  if(gamma <= 0 | gamma > 1)
    stop("p-value improvement factor gamma must be taken in ]0;1]")
  if(length(Tw) < wmax)
    stop("length of threshold vector Tw must be greater or equal to wmax")

  Raw.segments <- c()
  # ---------------------------------------------------------------------------.
  capt.time <- proc.time()[3]
  # ---------------------------------------------------------------------------.
  Raw.segments <- .C(
    "localopt_fisher",
    as.double(Yi),        # INPUT: sequence of measurements
    as.integer(n.probes), # INPUT: number of measurements
    as.character(name),   # INPUT: name used for the output file
    as.integer(wmax),     # INPUT: maximal scale to explore as number of measurements
    as.integer(wmin),	  # INPUT: minimal domain size as number of measurements
    as.double(gamma),	  # INPUT: p-value improvement factor
    as.double(Tw),        # INPUT: p-value threshold for each scale w
    Pw.min = double(n.probes),		 # OUTPUT: best p-value for each scale w
    n.segments = integer(1)  # OUTPUT: number of optimal domains
  )

  # ---------------------------------------------------------------------------.
  time.optimization <- proc.time()[3] - capt.time
  # ---------------------------------------------------------------------------.
  time.MRT.Analysis <- 0

  Raw.segments$file <- paste(name, "_RawSegments.txt", sep="")
  result <- list(
    n.segments      = 0,
    n.domains       = 0,
    n.maxresolution = 0,
    n.maxscale      = 0
  )
  if(length(Raw.segments)>0) {
    result <- list(
      time.optimization  = time.optimization,
      n.segments      = Raw.segments$n.segments,
      n.domains       = 0,
      n.maxresolution = 0,
      n.maxscale      = 0,
      file.segments   = Raw.segments$file,
      Pw.min          = Raw.segments$Pw.min
    )
    if(Raw.segments$n.segments>0) {
      # Multi Resolution Tree Analysis
      capt.time <- proc.time()[3]
      Domains <- .C(
        "MRT_Analysis",
        as.integer(n.probes), # INPUT: number of measurements
        as.character(name),   # INPUT: name used for the output file
        n.all      = integer(1), # OUTPUT: total number of simplified domains
        n.smallest = integer(1), # OUTPUT: number of Maximum Resolution Domains
        n.largest  = integer(1)  # OUTPUT: number of Maximum Scale Domains
      )
      time.MRT.Analysis <- proc.time()[3] - capt.time
      result <- list(
        time.optimization  = time.optimization,
        time.MRT.Analysis  = time.MRT.Analysis,
        n.segments         = Raw.segments$n.segments,
        n.domains          = Domains$n.all,
        n.maxresolution    = Domains$n.smallest,
        n.maxscale         = Domains$n.largest,
        file.segments      = Raw.segments$file,
        file.domains       = paste(name,"_Domains.txt",sep=""),
        file.maxresolution = paste(name,"_MaxResolutionDomains.txt",sep=""),
        file.maxscale      = paste(name,"_MaxScaleDomains.txt",sep=""),
        Pw.min            = Raw.segments$Pw.min
      )
    }
  }
  result
}
