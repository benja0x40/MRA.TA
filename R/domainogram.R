# =============================================================================.
# 
# -----------------------------------------------------------------------------.
calc.Qi <- function(x) {
  n <- length(x)
  r <- n-rank(x)+1
  q <- (r-0.5)/n
  q
}
# =============================================================================.
# 
# -----------------------------------------------------------------------------.
enrichmentScore <- function(x) {
  log(calc.Qi(x))
}
# =============================================================================.
# 
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
# 
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
# 
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
# 
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
# 
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
# 
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
