# =============================================================================.
#' Background bias estimation
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export backgroundBiasEstimation
#' @seealso
#'    \link{backgroundBiasCorrection},
#'    \link{normalizeArrayData}
# -----------------------------------------------------------------------------.
#' @description
#' Estimate a background bias within A, M values from a tiling array.
#'
##' @references
#' Leblanc B., Comet I., Bantignies F., and Cavalli G.,
#' Chromosome Conformation Capture on Chip (4C): data processing.
#' Book chapter in Polycomb Group Proteins: Methods and Protocols.
#' Lanzuolo C., Bodega B. editors, Methods in Molecular Biology (2016).
#' \link{http://dx.doi.org/10.1007/978-1-4939-6380-5_21}
# -----------------------------------------------------------------------------.
#' @param x
#' numeric vector of A values, e.g. the average of log2 of the specific and the
#' reference signals.
#'
#' @param y
#' numeric vector of M values, e.g. the log2 ratio of the specific over the
#' reference signal.
#'
#' @param AM.scale.compensation
#' \code{logical} determining if A and M values should be rescaled prior to
#' bias estimation. Scale compensation is required and activated by default.
#'
#' @param smoothness
#' numeric coefficient controlling the meanshift radius  in the rescaled (A,M)
#' plane.
#'
#' @param epsilon
#' convergence limit corresponding to the minimum displacement required per
#' meanshift iteration in the rescaled (A,M) plane.
#'
#' @param nsteps
#' minimum number of meanshift iterations.
#'
#' @param plots
#' \code{logical} used to activate control plots (default = F, no plots)
#'
#' @param xlab
#' label of the x axis for control plots.
#'
#' @param ylab
#' label of the y axis for control plots.
#'
# -----------------------------------------------------------------------------.
#' @return
#' numeric value of the estimated bias angle in radians
# -----------------------------------------------------------------------------.
#' @examples
#' # A simple test of the accuracy of bias estimations
#' ntst <- 10 # Increase to over 1000 for a reliable test of accuracy
#'
#' rotate <- function(x, y, theta) {
#'   u <- x * cos(theta) - y * sin(theta)
#'   v <- x * sin(theta) + y * cos(theta)
#'   list(x=u, y=v)
#' }
#'
#' tst <- c()
#' plots <- T
#' pb <- txtProgressBar(min=1, max=ntst, char="|", style=3)
#' for(i in 1:ntst) {
#'   n.bg <- sample(c(70000, 80000, 90000), 1) # Background level data
#'   n.sp <- 100000 - n.bg                     # Enriched level data
#'   a <- c(rnorm(n.bg, 10, 2),  rnorm(n.sp, 10, 1))
#'   m <- c(rnorm(n.bg, 0,  0.5), rnorm(n.sp,  1.5, 1))
#'   alpha <- runif(1, -pi/4, pi/4)            # Pick a random angle
#'   r <- rotate(a, m, theta = alpha)          # Simulate bias
#'   xtm <- Sys.time()
#'   theta <- backgroundBiasEstimation(        # Estimate bias
#'     r$x, r$y, plots=plots, AM.scale.compensation = F, smoothness = 0.07
#'   )
#'   xtm <- Sys.time() - xtm
#'   alpha <- 180/pi * alpha
#'   theta <- 180/pi * theta
#'   tst <- rbind(
#'     tst, c(n.bg=n.bg, n.sp=n.sp, time=xtm, simulated=alpha, estimated=theta)
#'   )
#'   setTxtProgressBar(pb, i)
#'   plots <- F
#' }
#' close(pb)
#' write.table(
#'   tst, "testBackgroundBiasEstimation.txt",
#'   quote = F, sep = '\t', row.names=F, col.names=T
#' )
#'
#' tst <- read.delim("testBackgroundBiasEstimation.txt", stringsAsFactors=F)
#' boxplot(tst$estimated - tst$simulated, range=0, ylab="estimated - simulated (degrees)")
#' legend("topright", paste("N =", ntst), bty="n")
#'
#' # Result with ntst=1000, delta <2.5 degrees in 95% of the tests:
#' # delta <- abs(tst$estimated - tst$simulated)
#' # q <- quantile(delta, probs=c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm=T)
# -----------------------------------------------------------------------------.
backgroundBiasEstimation <- function(x, y, AM.scale.compensation=T, smoothness=0.08, epsilon=0.001, nsteps=11, plots=F, xlab='A', ylab='M') {

  # Discard NA or infinite values
  discard <- (is.na(x) | is.na(y) | is.infinite(x) | is.infinite(y))
  ok.idx <- which(! discard)
  x <- x[ok.idx]
  y <- y[ok.idx]

  # Scale compensation (required for A and M values)
  if(AM.scale.compensation) {
    x <- x * sqrt(2)
    y <- y * sqrt(2)/2
  }
  DATA <- cbind(x,y)

  # Compute maximal amplitude of scale compensated (A,M) values
  maxamp <- sqrt((max(x)-min(x))^2 + (max(y)-min(y))^2)
  radius <- smoothness*maxamp

  # Pre-compute kd-tree structure of the (A,M) values for faster execution
  kdt <- kdtree(DATA, 1:nrow(DATA), maxamp/50)

  if(plots) {
    # MA-plot 1 (overlayed with a visualization of the background traversal)
    plot(
      x, y, xlab=xlab, ylab=ylab, pch='.', col=rgb(0,0,0,0.1)
    )
    # plot.kdtree.border(kdt)
    legend(
      "topleft",
      c("meanshift radius", "background traversal", "estimated bias"),
      fill=c(rgb(0,1,0), rgb(1,0.5,0), rgb(0.5,0,1)),
      bty='n'
    )
  }

  # Stage 1: find a good starting point (probe with minimal A value)
  i <- which.min(x)
  v0 <- c(x[i],y[i])
  for(i in 1:5) {
    # if(plots) points(v0[1], v0[2], pch='+', col=rgb(0,0,0,0.25))
    IW <- kdtree.getDN(DATA, v0, d=radius, node=kdt)
    n <- nrow(IW)
    if(length(n)==0) break
    v0[2] <- mean(DATA[IW[,1],2])
  }

  # Stage 2: background traversal using meanshift iterations
  bgwalk <- mean.shift(
    DATA, v0, d=radius, node=kdt,
    epsilon=epsilon, nsteps=nsteps, plots=plots
  )
  ctr   <- bgwalk$centroid
  alpha <- bgwalk$alpha

  # Stage 3: estimate bias slope as a trimmed mean of the traversal angles
  ra <- rank(alpha)
  a1 <- floor(length(alpha)/4)
  a2 <- floor(3*length(alpha)/4)
  theta <- mean(alpha[ra>=a1 & ra<=a2])

  if(plots) {
    # Superimpose visuzalization of estimated bias slope
    sbb <- line.ends(x, y, ctr, theta)
    segments(sbb$x1, sbb$y1, sbb$x2, sbb$y2, col=rgb(0.5,0,1), lwd=2)
    # Redraw background traversal
    moves <- bgwalk$moves
    nmvs  <- nrow(moves)
    last  <- which(ra>=a1 & ra<=a2)
    last  <- last[length(last)-1]
    segments(moves[(2:nmvs-1),1], moves[(2:nmvs-1),2], moves[2:nmvs,1], moves[2:nmvs,2], col=rgb(1,0.5,0),lwd=2)
    arrows(moves[last,1], moves[last,2], ctr[1], ctr[2], col=rgb(1,0.5,0),length=0.11,lwd=2)
    # Visualize center and radius of meanshift at starting and ending points
    # points(v0[1], v0[2], pch='+', col=rgb(1,0,0))
    # draw.circle(v0[1],v0[2],radius,nv=100,border=rgb(0,1,0),col=NA,lty=1,lwd=2)
    # points(ctr[1],ctr[2], pch='+', col=rgb(1,0,0))
    draw.circle(ctr[1],ctr[2],radius,nv=100,border=rgb(0,1,0),col=NA,lty=1,lwd=2)
    # Visualize the whole series of alpha angles and estimated bias slope
    plot(c(1,length(alpha)), c(-90,90), col="white", xlab="iteration", ylab="angle", main="Bias estimation")
    segments(0,theta*180/pi, length(alpha), theta*180/pi,col=rgb(0.5,0,1), lwd=2)
    points(alpha/pi*180,type='l', col=rgb(1,0.5,0), lwd=2)
    points((1:length(alpha))[ra>=a1 & ra<=a2],(alpha/pi*180)[ra>=a1 & ra<=a2],pch='+', col=rgb(0.5,0.5,0.5,0.75))

    text(length(alpha),theta/pi*180,labels=paste(round(theta/pi*180,1),"°",sep=""),col=rgb(0.5,0,1),pos=1)
  }

  theta
}
