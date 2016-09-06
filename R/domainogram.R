# =============================================================================.
#' Rank based scores
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export calc.Qi
#' @seealso
#'    \link{enrichmentScore},
#'    \link{domainogram},
#'    \link{segmentation}
# -----------------------------------------------------------------------------.
#' @description
#' Transform numeric values into rank-based scores.
#'
#' @references
#' de Wit E., Braunschweig U., Greil F., Bussemaker H.J., and van Steensel B.
#' Global chromatin domain organization of the Drosophila genome.
#' PLoS genetics 4: e1000045 (2008).
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/18369463}
# -----------------------------------------------------------------------------.
#' @param x vector of numeric values.
# -----------------------------------------------------------------------------.
#' @return
#' calc.Qi(x) = (rank(x) - 0.5) / n\cr where n is the length of x and ranking
#' is performed in decreasing order:
#' \itemize{
#'    \item rank(min(x)) = n
#'    \item rank(max(x)) = 1
#' }
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
calc.Qi <- function(x) {
  n <- length(x)
  r <- n-rank(x)+1
  q <- (r-0.5)/n
  q
}

# =============================================================================.
#' Log-transformed enrichment scores
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export enrichmentScore
#' @seealso
#'    \link{calc.Qi},
#'    \link{domainogram},
#'    \link{segmentation}
# -----------------------------------------------------------------------------.
#' @description
#' Compute log-transformed rank-based scores from numeric values.
# -----------------------------------------------------------------------------.
#' @param x vector of numeric values
# -----------------------------------------------------------------------------.
#' @return log(calc.Qi(x))
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
enrichmentScore <- function(x) {
  log(calc.Qi(x))
}

# =============================================================================.
#' Continuous coordinates for domainogram visualization (x axis)
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export visualizationCoordinates
#' @seealso
#'    \link{domainogram},
#'    \link{wSize2yAxis}
# -----------------------------------------------------------------------------.
#' @description
#' Compute coordinates for the domainogram visualization of a range of genomic
#' positions (x axis).
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param start
#' vector defining the start positions of genomic intervals.
#'
#' @param end
#' vector defining the end positions of genomic intervals.
# -----------------------------------------------------------------------------.
#' @return
#' \item{start }{start positions for domainogram visualization}
#' \item{end }{end positions for domainogram visualization}
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
visualizationCoordinates <- function(start, end) {
  x.s <- start
  x.e <- end
  idx <- 2:length(x.s)
  x.m <- (x.s[idx] + x.e[idx-1])/2
  x.e[(idx-1)] <- x.m
  x.s[idx] <- x.m
  list(start=x.s, end=x.e)
}

# =============================================================================.
#' Continuous coordinates for domainogram visualization (y axis)
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export wSize2yAxis
#' @seealso
#'    \link{domainogram},
#'    \link{visualizationCoordinates},
# -----------------------------------------------------------------------------.
#' @description
#' Compute coordinates for the domainogram visualization of a range of window
#' sizes (y axis).
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param wmax
#' integer indicating the maximum window size to be represented
#'
#' @param logscale
#' logical indicating if window sizes should be represented according to a log2
#' scale instead of a linear scale (default = F, linear)
# -----------------------------------------------------------------------------.
#' @return
#' \item{logscale}{value of the \code{logscale} parameter}
#' \item{ylab}{text to be used for the y axis legend}
#' \item{y.t}{top coordinate for each window size}
#' \item{y.b}{bottom coordinate for each window size}
#' \item{ylim}{range of the visualization coordinates (y axis)}
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
wSize2yAxis <- function(wmax, logscale=FALSE) {
  res <- list()
  res$logscale <- logscale
  res$ylab <- 'window size'
  res$y.b <- 1:wmax-1
  res$y.t <- 1:wmax
  if(logscale) {
    res$ylab <- 'log2(window size)'
    res$y.b <- log2(1:wmax)
    res$y.t <- log2(2:(wmax+1))
  }
  res$ylim <- c(min(res$y.b), max(res$y.t))
  res
}
# =============================================================================.
#' Create a new multi-resolution plot
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export
#' @seealso
#'    \link{domainogram},
#'    \link{plotOptimalSegments}
# -----------------------------------------------------------------------------.
#' @description
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param x.s
#'
#' @param x.e
#'
#' @param w2y
#'
#' @param xlim
#'
#' @param wlim
#'
#' @param xlab
#'
# -----------------------------------------------------------------------------.
#' @return NULL
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
new.domain.plot <- function(x.s, x.e, w2y, xlim=NULL, wlim=NULL, xlab='') {
  if(is.null(xlim)) {
    xlim <- c(min(x.s), max(x.e))
  }
  if(is.null(wlim)) {
    ylim <- w2y$ylim
  } else {
    ylim <- c(w2y$y.b[wlim[1]], w2y$y.b[wlim[2]])
  }
  message("[new.domain.plot]")
  message("xlim = ", paste(xlim, collapse=" "))
  message("ylim = ", paste(ylim, collapse=" "))
  message("wlim = ", paste(wlim, collapse=" "))
  message("logscale = ", paste(w2y$logscale, collapse=" "))
  plot(0, xlim=xlim, ylim=ylim, type='n', xlab=xlab, ylab=w2y$ylab, col='white', axes=FALSE, xaxs='i',  yaxs='i')
  axis(1)
  segments(xlim[1], ylim[1], xlim[2], ylim[1], lwd=1)
  axis(2)
  segments(xlim[1], ylim[1], xlim[1], ylim[2], lwd=1)
}

# =============================================================================.
#' Plot multi-resolution segmentation results
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export plotOptimalSegments
#' @seealso
#'    \link{visualizationCoordinates},
#'    \link{wSize2yAxis},
#'    \link{domainogram},
#'    \link{segmentation},
#'    \link{plotDomains}
# -----------------------------------------------------------------------------.
#' @description
#' Plot optimal segments resulting from the \link{segmentation} function.
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param optimums
#' data.frame of segmentation results generated by the \link{segmentation}
#' function.
#'
#' @param x.s
#' vector of start coordinates (x axis) defining the represented genomic
#' intervals.
#'
#' @param x.e
#' vector of end coordinates (x axis) defining the represented genomic
#' intervals.
#'
#' @param w2y
#' coordinate mapping for window sizes (y axis) produced by the
#' \link{wSize2yAxis} function.
#'
#' @param xlim
#' range of the represented genomic region, which must be indicated as
#' c(xmin, xmax).
#'
#' @param wlim
#' range of represented window sizes, which must be indicated as c(wmin, wmax).
#'
#' @param xlab
#' legend of the x axis.
#'
#' @param add
#' logical, when set to TRUE the plot should overlay an existing multi-resolution
#' plot.
#'
#' @param type
#' either "windows" (default) or "middles".
#'
#' @param \dots
#' optional parameters forwarded to the \code{rect} function.
# -----------------------------------------------------------------------------.
#' @return NULL
# -----------------------------------------------------------------------------.
#' @examples
#' # Simulate enrichment signal
#' n <- 2000
#' Mi <- rep(0, n)
#' Mi <- Mi + dnorm(1:n, 2.5*n/20,  n/40) + dnorm(1:n, 4*n/20,  50)
#' Mi <- Mi + 4 * dnorm(1:n, 5*n/10, n/10)
#' Mi <- Mi + dnorm(1:n, 16*n/20,  n/40) + dnorm(1:n, 17.5*n/20,  50)
#' Mi <- (Mi/max(Mi))^4 + rnorm(n)/4
#' # Compute enrichment scores
#' Yi <- enrichmentScore(Mi)
#' # Multi-resolution segmentation
#' seg.c <- segmentation(Yi, name="MRA_demo", wmin=20)
#' # Load segmentation results
#' opts <- read.delim(seg.c$file.segments, stringsAsFactors=F, skip=1)
#' doms <- read.delim(seg.c$file.domains, stringsAsFactors=F)
#' doms.mr <- read.delim(seg.c$file.maxresolution, stringsAsFactors=F)
#' doms.ms <- read.delim(seg.c$file.maxscale, stringsAsFactors=F)
#' # Visualization coordinates
#' x.s <- 1:n - 0.5
#' x.e <- 1:n + 0.5
#' w2y <- wSize2yAxis(n, logscale=T)
#' layout(matrix(1:2, 2, 1), heights=c(3,1)/4)
#' par(mar=c(3, 4, 1, 2)) # default bottom, left, top, right = c(5, 4, 4, 2)
#' # Plot domainogram
#' domainogram(Yi, x.s, x.e, w2y)
#' plot(Mi, type='l')
#' # Visualize segmentation results
#' plotOptimalSegments(opts, x.s, x.e, w2y, col="black")
#' plot(Mi, type='l')
#' # Visualize multi-resolution domains
#' plotDomains(doms, x.s, x.e, w2y, col=rgb(0,0,0,0.5), border=rgb(0,0,0,0))
#' # Visualize max. resolution and max. scale domains
#' plotDomains(doms.mr, x.s, x.e, w2y, col=rgb(0,1,0,0.5), border=rgb(0,0,0,0), lwd=2, lty=1, add=T)
#' plotDomains(doms.ms, x.s, x.e, w2y, col=rgb(1,0,0,0.5), border=rgb(0,0,0,0), lwd=2, lty=1, add=T)
#' legend("topright", c("Max. resolution", "Max. scale", "Both"), fill=c("green", "red", "chocolate"), bty='n')
#' plot(Mi, type='l')
# -----------------------------------------------------------------------------.
plotOptimalSegments <- function(optimums, x.s, x.e, w2y, xlim=NULL, wlim=NULL, xlab='', add=F, type=c("windows", "middles"), ...) {

  type <- type[1]

  if(! add) { # Make new plot
    new.domain.plot(x.s, x.e, w2y, xlim, wlim, xlab)
  }

  if(nrow(optimums)>0) {
    # Skip optimums that are entirely invisible         # *** WARNING/TODO *** #
    if(type=="windows") {
      rect(x.s[optimums$i-optimums$w+1], w2y$y.b[optimums$w], x.e[optimums$i], w2y$y.t[optimums$w], ...)
    }
    if(type=="middles") {
      i <- ceiling(optimums$i-optimums$w/2)
      rect(x.s[i], w2y$y.b[optimums$w], x.e[i], w2y$y.t[optimums$w], ...)
    }
  }
}

# =============================================================================.
#' Plot multi-resolution segmentation results
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export plotDomains
#' @seealso
#'    \link{visualizationCoordinates},
#'    \link{wSize2yAxis},
#'    \link{domainogram},
#'    \link{segmentation},
#'    \link{plotOptimalSegments}
# -----------------------------------------------------------------------------.
#' @description
#' Plot multi-resolution domains resulting from the \link{segmentation}
#' function.
#'
#' @details
# -----------------------------------------------------------------------------.
#' @param domains
#' data.frame of segmentation results generated by the \link{segmentation}
#' function.
#'
#' @param x.s
#' vector of start coordinates (x axis) defining the represented genomic
#' intervals.
#'
#' @param x.e
#' vector of end coordinates (x axis) defining the represented genomic
#' intervals.
#'
#' @param w2y
#' coordinate mapping for window sizes (y axis) produced by the
#' \link{wSize2yAxis} function.
#'
#' @param xlim
#' range of the represented genomic region, which must be indicated as
#' c(xmin, xmax).
#'
#' @param wlim
#' range of represented window sizes, which must be indicated as c(wmin, wmax).
#'
#' @param xlab
#' legend of the x axis.
#'
#' @param add
#' logical, when set to TRUE the plot should overlay an existing multi-resolution
#' plot.
#'
#' @param \dots
#' optional parameters forwarded to the \code{rect} function.
# -----------------------------------------------------------------------------.
#' @return NULL
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
#' plot(Mi, type='l')
#'
#' # Visualize segmentation results
#' plotOptimalSegments(opts, x.s, x.e, w2y, col="black")
#' plot(Mi, type='l')
#'
#' # Visualize multi-resolution domains
#' plotDomains(doms, x.s, x.e, w2y, col=rgb(0,0,0,0.5), border=rgb(0,0,0,0))
#' # Visualize max. resolution and max. scale domains
#' plotDomains(doms.mr, x.s, x.e, w2y, col=rgb(0,1,0,0.5), border=rgb(0,0,0,0), lwd=2, lty=1, add=T)
#' plotDomains(doms.ms, x.s, x.e, w2y, col=rgb(1,0,0,0.5), border=rgb(0,0,0,0), lwd=2, lty=1, add=T)
#' legend("topright", c("Max. resolution", "Max. scale", "Both"), fill=c("green", "red", "chocolate"), bty='n')
#' plot(Mi, type='l')
# -----------------------------------------------------------------------------.
plotDomains <- function(domains, x.s, x.e, w2y, xlim=NULL, wlim=NULL, xlab='', add=F, ...) {

  if(! add) { # Make new plot
    new.domain.plot(x.s, x.e, w2y, xlim, wlim, xlab)
  }

  if(nrow(domains)>0) {
    # Skip domains that are entirely invisible          # *** WARNING/TODO *** #
    rect(x.s[domains$start], w2y$y.b[domains$wmin], x.e[domains$end], w2y$y.t[domains$wmax], ...)
  }
}

# =============================================================================.
#' Domainogram visualization
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export domainogram
#' @seealso
#'    \link{calc.Qi},
#'    \link{enrichmentScore},
#'    \link{visualizationCoordinates},
#'    \link{wSize2yAxis},
#'    \link{makeColors},
#'    \link{segmentation},
#'    \link{plotOptimalSegments},
#'    \link{plotDomains}
# -----------------------------------------------------------------------------.
#' @description
#' Plot a domainogram from a sequence of log-transformed rank-based scores.
#' This function allows to use a linear or log2 scale for window sizes (y
#' axis).
#'
#' @details
#' The domainogram is a multi-resolution representation of scan statistics
#' introduced by de Wit et al., 2008 (see references below). It is based on
#' R.A. Fisher's method for the combination of p-values applied to
#' non-parametric scores.
#'
#' @references
#' de Wit E., Braunschweig U., Greil F., Bussemaker H.J., and van Steensel B.
#' Global chromatin domain organization of the Drosophila genome.
#' PLoS genetics 4: e1000045 (2008).
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/18369463}
# -----------------------------------------------------------------------------.
#' @param Yi
#' vector of statistical scores produced by the \link{enrichmentScore} function
#'
#' @param x.s
#' vector of start coordinates (x axis) defining the represented genomic
#' intervals.
#'
#' @param x.e
#' vector of end coordinates (x axis) defining the represented genomic
#' intervals.
#'
#' @param w2y
#' coordinate mapping for window sizes (y axis) produced by the
#' \link{wSize2yAxis} function.
#'
#' @param xlim
#' range of the represented genomic region, which must be indicated as
#' c(xmin, xmax).
#'
#' @param wlim
#' range of represented window sizes, which must be indicated as c(wmin, wmax).
#'
#' @param colfun
#' coloring function. The function \link{makeColors} is used by default.
#'
#' @param xlab
#' legend of the x axis.
#'
#' @param add
#' logical, when set to TRUE the plot should overlay an existing multi-resolution
#' plot.
#'
#' @param logskip
#' numeric coefficient for skipping the representation of window sizes when
#' using the log2 scale on the y axis.
#' A value of 1.0 forces to represent all window sizes.
#' A value of 2 would double the window sizes represented at each iteration of
#' the rendering.
#'
#' @param \dots
#' optional parameters forwarded to the \code{colfun} function.
#' When using the default \link{makeColors} function, these parameters are
#' \code{thresholds}, \code{colors}, \code{background} and \code{overflow}
#' (see documentation for details).
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
#' # Simulate enrichment signal
#' n <- 500
#' Mi <- cos(2*pi*2.5*(-n:n)/n) + rnorm(2*n+1)
#' n <- length(Mi)
#'
#' # Compute enrichment scores
#' Yi <- enrichmentScore(Mi)
#'
#' # Visualization coordinates
#' x.s <- 1:n - 0.5
#' x.e <- 1:n + 0.5
#' w2y <- wSize2yAxis(n)
#'
#' layout(matrix(1:2, 2, 1), heights=c(3,1)/4)
#' par(mar=c(3, 4, 1, 2)) # default bottom, left, top, right = c(5, 4, 4, 2)
#'
#' # Plot domainogram
#' domainogram(Yi, x.s, x.e, w2y)
#' plot(Mi, type='l')
# -----------------------------------------------------------------------------.
domainogram <- function (Yi, x.s, x.e, w2y, xlim=NULL, wlim=NULL, colfun=makeColors, xlab='', add=F, logskip=1.01, ...) {

  if(! add) { # Make new plot
    new.domain.plot(x.s, x.e, w2y, xlim, wlim, xlab)
  }
  if(is.null(xlim)) {
    xlim <- c(min(x.s), max(x.e))
  }
  # Count number of measurements and define the range of window sizes
  n.probes <- length(Yi)
  if(is.null(wlim)) {
    wlim <- c(1, n.probes)
  }
  wmin <- wlim[1]
  wmax <- wlim[2]

  # Visible windows at w = wmin
  idx <- which(x.e>=xlim[1] & x.s<=xlim[2])
  x1 <- x.s[idx]
  x2 <- x.e[idx]

  y1 <- rep(w2y$y.b[wmin], length(idx))
  y2 <- rep(w2y$y.t[wmin+1], length(idx))
  rect(x1, y1, x2, y2, col=colfun(Yi[idx], ...), border=rgb(0,0,0,0))

  # Initialization
  S1 <- Yi[1]
  winc <- 1
  w <- 2

  # Domainogram (visualization of multi resolution statistics)
  pb <- txtProgressBar(min = 1, max = wmax-1, style = 3, char="|")
  while(w <= wmax) {

    # Compute windowed p-values using R.A. Fisher's combine method
    S1 <- S1 + sum(Yi[(w-winc+1):w])
    Pw <- c(S1, diff(Yi, w))
    Pw <- cumsum(Pw)
    Pw <- pchisq(-2*Pw, df=2*w, lower.tail=FALSE, log.p=TRUE)
    wn <- length(Pw)

    # Plot domainogram
    if(w >= wmin) {
      if(w%%2!=0) {
        x1 <- x.s[floor(w/2)+(1:wn)]
        x2 <- x.e[floor(w/2)+(1:wn)]
      }
      else {
        x1 <- (x.s[floor(w/2)+(1:wn)-1]+x.e[floor(w/2)+(1:wn)-1])/2
        x2 <- (x.s[floor(w/2)+(1:wn)]+x.e[floor(w/2)+(1:wn)])/2
      }
      idx <- which(x2>=xlim[1] & x1<=xlim[2])
      x1 <- x1[idx]
      x2 <- x2[idx]
      y1 <- rep(w2y$y.b[w - winc + 1], length(idx))
      y2 <- rep(w2y$y.t[w], length(idx))
      rect(x1, y1, x2, y2, col=colfun(Pw[idx], ...), border=rgb(0,0,0,0))
    }

    # Increase window size
    if(w2y$logscale) {
      winc <- max(1, round(1.01*w)-w)
    }
    w <- w + winc

    setTxtProgressBar(pb, w)
  }
  close(pb)
}
