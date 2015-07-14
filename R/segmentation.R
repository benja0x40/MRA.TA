# =============================================================================.
# 
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
