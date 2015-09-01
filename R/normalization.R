# FUNCTIONS ####################################################################

# =============================================================================.
# FUNCTION line.ends
# DETAIL ----------------------------------------------------------------------.
# Utility function for visualization of the (A,M) bias slope.
# (see normalizeArrayData function below)
# -----------------------------------------------------------------------------.
# USAGE NOT DOCUMENTED
# -----------------------------------------------------------------------------.
line.ends <- function(x, y, ctr, alpha) {
  x1a <- min(x, na.rm=TRUE)
  x2a <- max(x, na.rm=TRUE)
  y1a <- ctr[2] - sin(alpha) * (ctr[1] - x1a) / cos(alpha)
  y2a <- ctr[2] + sin(alpha) * (x2a - ctr[1]) / cos(alpha)

  y1b <- min(y, na.rm=TRUE)
  y2b <- max(y, na.rm=TRUE)
  x1b <- ctr[1] - cos(alpha) * (ctr[2] - y1b) / sin(alpha)
  x2b <- ctr[1] + cos(alpha) * (y2b - ctr[2]) / sin(alpha)

  M <- c()
  if(x1b>=x1a & x1b<=x2a) M <- rbind(M, c(x1b, y1b))
  if(x2b>=x1a & x2b<=x2a) M <- rbind(M, c(x2b, y2b))
  if(y1a>=y1b & y1a<=y2b) M <- rbind(M, c(x1a, y1a))
  if(y2a>=y1b & y2a<=y2b) M <- rbind(M, c(x2a, y2a))

  if(M[1,1]<M[2,1])
    return(list(x1=M[1,1], y1=M[1,2], x2=M[2,1], y2=M[2,2]))
  else
    return(list(x1=M[2,1], y1=M[2,2], x2=M[1,1], y2=M[1,2]))
}
# =============================================================================.
# FUNCTION plot.kdtree.border
# DETAIL ----------------------------------------------------------------------.
# Utility function for visualizing the kd-tree structure of (A,M) values which
# is used to accelerate the background traversal algorithm.
# (see normalizeArrayData function below)
# -----------------------------------------------------------------------------.
# USAGE NOT DOCUMENTED
# -----------------------------------------------------------------------------.
plot.kdtree.border <- function(node) {
  if(node$sg.tag) {
    plot.kdtree.border(node$sg.inf)
    plot.kdtree.border(node$sg.sup)
    rect(node$sg.min[1], node$sg.min[2], node$sg.max[1], node$sg.max[2], col=rgb(0,0,1,0.02), border=rgb(0,0,0,0.05))
  }
  else {
    rect(node$sg.min[1], node$sg.min[2], node$sg.max[1], node$sg.max[2], col=rgb(0,0,1,0.02), border=rgb(0,0,0,0.05))
  }
}
# =============================================================================.
# FUNCTION kdtree
# DETAIL ----------------------------------------------------------------------.
# Utility function generating a kd-tree structure of (A,M) values for faster
# execution of the background traversal algorithm.
# (see normalizeArrayData function below)
# -----------------------------------------------------------------------------.
# USAGE NOT DOCUMENTED
# -----------------------------------------------------------------------------.
kdtree <- function(M, i, d) {

  n <- length(i)
  node <- list()

  node$sg.tag <- FALSE
  node$sg.nbr <- n

  if(n==1) {
    node$sg.min <- M[i,]
    node$sg.max <- M[i,]
    node$sg.avg <- M[i,]
    node$sg.ids <- i
    return(node)
  }

  node$sg.min <- apply(M[i,], MARGIN=2, FUN=min)
  node$sg.max <- apply(M[i,], MARGIN=2, FUN=max)
  node$sg.avg <- apply(M[i,], MARGIN=2, FUN=mean)

  s.d <- abs(node$sg.max - node$sg.min)
  s.c <- as.integer(which.max(s.d))

  if( s.d[s.c] < d) {
    node$sg.ids <- i
  }
  else {
    s.m <- median(M[i, s.c])

    node$sg.tag <- TRUE
    node$sg.cut <- s.c
    if(s.m - min(M[i, s.c]) > max(M[i, s.c]) - s.m ) {
      node$sg.inf <- kdtree(M, i[ M[i, s.c] <  s.m ], d)
      node$sg.sup <- kdtree(M, i[ M[i, s.c] >= s.m ], d)
    }
    else {
      node$sg.inf <- kdtree(M, i[ M[i, s.c] <= s.m ], d)
      node$sg.sup <- kdtree(M, i[ M[i, s.c] >  s.m ], d)
    }
  }
  return(node)
}
# =============================================================================.
# FUNCTION kdtree.getDN
# DETAIL ----------------------------------------------------------------------.
# Utility function for fast neighbor selection during the background traversal
# algorithm, based on a pre-computed kd-tree structure of (A,M) values.
# (see normalizeArrayData function below)
# -----------------------------------------------------------------------------.
# USAGE NOT DOCUMENTED
# -----------------------------------------------------------------------------.
kdtree.getDN <- function(M, v, d, node) {

  n <- node$sg.nbr
  m <- length(v)

  if(node$sg.tag) {
    a <- c()
    b <- c()
    s.c <- node$sg.cut

    if((node$sg.inf$sg.max[s.c] + d) >= v[s.c]) {
      a <- kdtree.getDN(M, v, d, node$sg.inf)
    }
    if((node$sg.sup$sg.min[s.c] - d) <= v[s.c]) {
      b <- kdtree.getDN(M, v, d, node$sg.sup)
    }
    return(rbind(a, b))
  }

  i <- node$sg.ids
  A <- M[i,] - matrix(v, n, m, byrow=TRUE)
  A <- sqrt(apply(A*A, MARGIN=1, FUN=sum))
  match <- A <= d
  return(cbind(i[ match ], A[ match ]))
}
# =============================================================================.
# FUNCTION mean.shift
# DETAIL ----------------------------------------------------------------------.
# Core function of the background traversal algorithm based on a pre-computed
# kd-tree structure of (A,M) values.
# (see normalizeArrayData function below)
# -----------------------------------------------------------------------------.
# USAGE NOT DOCUMENTED
# -----------------------------------------------------------------------------.
mean.shift <- function(M, v, d, node, epsilon=0.001, col=rgb(1,0.5,0,1), nsteps=0, plots=FALSE) {

  i <- rep(FALSE, nrow(M))
  m <- length(v)
  moves <- v
  alpha <- c()

  while(TRUE) {

    IW <- kdtree.getDN(M, v, d, node)

    n <- nrow(IW)

    if(length(n)==0)
      return(NULL)

    i[IW[,1]] <- TRUE

    A <- M[IW[,1],] - matrix(v, n, m, byrow=TRUE)
    A <- apply(A, MARGIN=2, FUN=mean)

    if(plots) segments(v[1], v[2], (v+A)[1], (v+A)[2], col=col,lwd=2)

    v <- v + A
    moves <- rbind(moves, v)

    a <- asin(A[2]/sqrt(A[1]*A[1]+A[2]*A[2]))
    alpha <- c(alpha, a)

    if(length(alpha)>nsteps & sqrt(sum(A*A)) < epsilon) break
  }

  list(i=i, centroid=v, alpha=alpha, moves=moves)
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
# x = A value
# y = M value
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
# =============================================================================.
# Utility function for planar rotation (2D) centered at origin (0,0)
# -----------------------------------------------------------------------------.
# USAGE NOT DOCUMENTED
# -----------------------------------------------------------------------------.
rotate <- function(x, y, theta) {
  u <- x * cos(theta) - y * sin(theta)
  v <- x * sin(theta) + y * cos(theta)
  list(x=u, y=v)
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
backgroundBiasCorrection <- function(x, y, theta, AM.scale.compensation=T) {
  # Scale compensation (required for A and M values)
  if(AM.scale.compensation) {
    x <- x * sqrt(2)
    y <- y * sqrt(2)/2
  }
  v <- rotate(x, y, -theta)
  v
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
# *** WARNING/TODO *** #  Cleanup the returned A value (and list names?)
lowessCorrection <- function(x, y, lowess.f=0.3, lowess.mad=0, lowess.iter=5, plots=F) {
  A <- x
  M <- y

  if(lowess.mad>0) {
    j <- which(abs(y-median(y))<lowess.mad*mad(y))
    normalization <- lowess(x[j], y[j], f=lowess.f, iter=lowess.iter)
    o <- order(A)
    yn <- approx(normalization$x, normalization$y, xout=A[o], rule=2, ties='ordered')$y
    y <- M - yn[order(o)]
    # x <- A + yn[order(o)]/2
  }
  else {
    j <- 1:length(x)
    normalization <- lowess(x, y, f=lowess.f, iter=lowess.iter)
    o <- order(A)
    yn <- normalization$y
    y <- M - yn[order(o)]
    # x <- A + yn[order(o)]/2
  }
  if(plots) {
    title <- paste("Lowess (f =", lowess.f)
    if(lowess.mad==0) title <- paste(title,")",sep="")
    else title <- paste(title, ", ", lowess.mad," x mad)", sep="")
    # MA-plot 4 (overlayed with a visualization of the lowess)
    plot(A, M ,pch='.', col=rgb(0,0,0,0.1), main=title)
    if(lowess.mad>0)
      points(A[j], M[j] ,pch='.', col=rgb(0,0.25,0.75,0.05))
    lines(A[o], yn, col=rgb(0.5,0,1,1))
    # MA-plot 5 (overlayed with a visualization of the lowess)
    plot(x, y ,pch='.', col=rgb(0,0,0,0.1), main="Lowess corrected", xlab='A', ylab='M')
    lines(c(min(x),max(x)),c(0,0),col=rgb(0.3,0.3,0.3,0.7),lwd=1.5)
  }
  list(x=x, y=y)
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
normalizeArrayData <- function(A, M, smoothness=0.08, epsilon=0.01, nsteps=11, name="Test", plots=TRUE, lowess=T, lowess.f=0.2, lowess.mad=0, lowess.iter=5) {

  # Discard NA or infinite values
  discard <- (is.na(A) | is.na(M) | is.infinite(A) | is.infinite(M))
  ok.idx <- which(! discard)
  if(length(ok.idx)!=length(A)) {
    message("Probes associated with at least one NA in {A,M} values will not be normalized")
  }
  x <- A[ok.idx]
  y <- M[ok.idx]

  # Graphic setup
  if(plots) {
    npix <- 500
    if(lowess) {
      png(paste(name, "NormalizationReport.png",sep="_"), width=3*npix ,height=2*npix)
      layout(matrix(1:6,2,3,byrow=T))
    }
    else {
      png(paste(name, "NormalizationReport.png",sep="_"), width=2*npix ,height=2*npix)
      layout(matrix(1:4,2,2,byrow=T))
    }
    par(cex=1.4)
  }

  # Background bias estimation and correction
  theta <- backgroundBiasEstimation(
    x, y, AM.scale.compensation=T, smoothness=smoothness, epsilon=epsilon, nsteps=nsteps, plots=plots, xlab='A', ylab='M'
  )
  v <- backgroundBiasCorrection(x, y, theta, AM.scale.compensation=T)
  x <- v$x
  y <- v$y

  # MA-plot 3 (bias corrected A and M values)
  if(plots) plot(x, y ,pch='.', xlab='A', ylab='M', col=rgb(0,0,0,0.1), main="Bias corrected")

  # Lowess correction
  if(lowess) {
    v <- lowessCorrection(x, y, lowess.f=lowess.f, lowess.mad=lowess.mad, lowess.iter=lowess.iter, plots=plots)
    x <- v$x
    y <- v$y
  }

  if(plots) {
    dev.off()
  }

  A[ok.idx] <- x
  M[ok.idx] <- y

  return(list(A=A, M=M, bias=theta, ignored=which(discard)))
}
# =============================================================================.
# Filter microarray probes according to simple quality scores
# -----------------------------------------------------------------------------.
# x = raw Ctr signal intensity
# y = raw 4C signal intensity
# -----------------------------------------------------------------------------.
enrichmentQuality <- function(x, y, q.best=0.015, q.worst=5E-3, plots=F) {

  n.probes <- length(x)

  A <- (log2(y) + log2(x))/2
  M <-  log2(y) - log2(x)

  qx <- rank(x)/n.probes
  qy <- rank(y)/n.probes
  qa <- rank(A)/n.probes
  qm <- rank(M)/n.probes
  # qs <- 1 - rank(abs(mean(range(A)) - A))/n.probes  # Not used

  f <- 1/2 # f = 0 => qb = qm
  qb <- (rank(qy - qx)/n.probes)^f * qm

  chk.best <- rank(qb)/n.probes > 1 - q.best

  qw <- qx * qa
  chk.worst <- rank(qw)/n.probes > 1 - q.worst
  if(plots) {
    clr <- rgb(1-qb, 1-qw, 0, 0.5)
    plot(x, y, pch='.', col=clr, main="Quality score (x,y)", log='xy')
    abline(0, 1, col="grey", lwd=2)

    clr <- rgb(chk.worst, chk.best, 0, ifelse(chk.worst | chk.best, 0.5, 0.1))
    plot(x, y, pch='.', col=clr, log='xy', main="Extrema (x,y)")
    abline(0, 1, col="grey", lwd=2)

    clr <- rgb(1-qb, 1-qw, 0, 0.5)
    plot(A, M, pch='.', col=clr, main="Quality score (A,M)")
    abline(h=0, col="grey", lwd=2)
    legend("topleft", "good", fill="green", bty='n')
    legend("topright", "poor", fill="red", bty='n')

    clr <- rgb(chk.worst, chk.best, 0, ifelse(chk.worst | chk.best, 0.5, 0.1))
    plot(A, M, pch='.', col=clr, main="Extrema (A,M)")
    abline(h=0, col="grey", lwd=2)
    legend("topleft", "best", fill="green", bty='n')
    legend("topright", "problematic", fill="red", bty='n')
  }

  list(is.best=chk.best, is.worst=chk.worst, q.best=q.best, q.worst=q.worst, b.score=qb, w.score=qw)
}

# =============================================================================.
#
# -----------------------------------------------------------------------------.
plotProbeDistanceControls <- function(x, rnk, dis, dlim=c(-2000, 2000), QF=NULL, n.q=100, ylab="4C") {
  n.c <- (n.q%/%2) - 1
  a <- 0:n.c/n.c
  b <- n.c:0/n.c
  clrs <- c(
    rgb(0.5*a, 0.5*a, 1 - 0.5*a),
    rgb(0.5+0.5*a, 0.5*b^(1/10), 0.5*b)
  )
  rank.ctrl <- tapply(x, INDEX=rnk, FUN=quantile, probs=0:n.q/n.q)
  rank.ctrl <- matrix(unlist(rank.ctrl), max(rnk, na.rm=T), n.q+1, byrow=T)
  rmax <- max(which(table(rnk, useNA="no")>=3*n.q))
  rmax <- min(rmax, nrow(rank.ctrl))
  plot(0, type='n', xlim=c(1, rmax), ylim=range(x), xaxs='i', xlab='probe rank', ylab=ylab, main = "percentiles")
  grid(nx = rmax, ny=0)
  abline(h=0)
  for(i in 1:ncol(rank.ctrl)) {
    lines(1:nrow(rank.ctrl), rank.ctrl[,i], col=clrs[i])
  }
  legend("topright", c("highest","median","lowest"), fill=c("red", "grey", "blue"), bty='n')
  #boxplot(x ~ rnk, range=1, outline=T, outpch=20, outcex=0.5, outcol=rgb(0,0,0,0.1), xlab='probe rank', ylab=ylab, xlim=c(0, nrow(rank.ctrl)+1), xaxs='i')
  #abline(h=0)
  clr <- rgb(0,0,0,0.1)
  pch='.' # 46
  if(! is.null(QF)) {
    clr <- rgb(QF$is.worst, QF$is.best, 0, ifelse(QF$is.worst | QF$is.best, 0.5, 0.1))
    pch = ifelse(QF$is.worst | QF$is.best, 20, 46)
  }
  plot(0, type='n', xlim=dlim, ylim=range(x, na.rm=T), xlab='probe distance', ylab=ylab)
  abline(h=0, col='lightgrey', lwd=2)
  abline(v=0, col=rgb(0.5,0,1), lwd=2)
  points(dis, x, pch=pch, col=clr, cex=0.5)
  if(! is.null(QF)) {
    legend("topleft", "best", fill="green", bty='n')
    legend("topright", "problematic", fill="red", bty='n')
  }
}

# =============================================================================.
#
# -----------------------------------------------------------------------------.
plotSelectedProbes <- function(a, m, dis, sel, dlim=c(-2000, 2000), ylab="4C") {
  clr <- rgb(0,0,0,0.1)
  pch=20
  clr <- rgb((1-sel)/1.1, sel/1.1, 0, 0.1)
  plot(0, type='n', xlim=dlim, ylim=range(m, na.rm=T), xlab='probe distance', ylab=ylab)
  abline(h=0, col='lightgrey', lwd=2)
  abline(v=0, col=rgb(0.5,0,1), lwd=2)
  # points(dis, m, pch=pch, col=rgb(0,0,0,0.1), cex=0.5)
  points(dis, m, pch=pch, col=clr, cex=0.5)
  legend("topleft", "selected", fill="green", bty='n')
  legend("topright", "rejected", fill="red", bty='n')
}

