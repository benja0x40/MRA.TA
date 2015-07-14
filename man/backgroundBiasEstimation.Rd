\name{backgroundBiasEstimation}
\alias{backgroundBiasEstimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Background bias estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Estimate a background bias between a specific and a reference signal.
}
\usage{
backgroundBiasEstimation(x, y, AM.scale.compensation = T, smoothness = 0.08, epsilon = 0.001, nsteps = 11, plots = F, xlab = "A", ylab = "M")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  %%     ~~Describe \code{x} here~~
  average log2 of the specific and the reference signal (A value)
}
  \item{y}{
  %%     ~~Describe \code{y} here~~
  log2 ratio of the specific over the reference signal (M value)
}
  \item{AM.scale.compensation}{
  %%     ~~Describe \code{AM.scale.compensation} here~~
}
  \item{smoothness}{
  %%     ~~Describe \code{smoothness} here~~
}
  \item{epsilon}{
  %%     ~~Describe \code{epsilon} here~~
}
  \item{nsteps}{
  %%     ~~Describe \code{nsteps} here~~
}
  \item{plots}{
  %%     ~~Describe \code{plots} here~~
}
  \item{xlab}{
  %%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
  %%     ~~Describe \code{ylab} here~~
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{backgroundBiasCorrection}},
  \code{\link{normalizeArrayData}}
}
\examples{
  ntst <- 10 # Increase to over 1000 for a reliable test of accuracy
  
  rotate <- function(x, y, theta) {
    u <- x * cos(theta) - y * sin(theta)
    v <- x * sin(theta) + y * cos(theta)
    list(x=u, y=v)
  }
  
  tst <- c()
  plots <- T
  pb <- txtProgressBar(min=1, max=ntst, char="|", style=3)
  for(i in 1:ntst) {
    n.bg <- sample(c(70000, 80000, 90000), 1) # Background level data
    n.sp <- 100000 - n.bg                     # Enriched level data
    a <- c(rnorm(n.bg, 10, 2),  rnorm(n.sp, 10, 1))
    m <- c(rnorm(n.bg, 0,  0.5), rnorm(n.sp,  1.5, 1))
    alpha <- runif(1, -pi/4, pi/4)            # Pick a random angle
    r <- rotate(a, m, theta = alpha)          # Simulate bias
    xtm <- Sys.time()
    theta <- backgroundBiasEstimation(        # Estimate bias
      r$x, r$y, plots=plots, AM.scale.compensation = F, smoothness = 0.07
    )
    xtm <- Sys.time() - xtm 
    alpha <- 180/pi * alpha
    theta <- 180/pi * theta
    tst <- rbind(
      tst, c(n.bg=n.bg, n.sp=n.sp, time=xtm, simulated=alpha, estimated=theta)
    )
    setTxtProgressBar(pb, i)
    plots <- F
  }
  close(pb)
  write.table(
    tst, "testBackgroundBiasEstimation.txt",
    quote = F, sep = '\t', row.names=F, col.names=T
  )
  
  tst <- read.delim("testBackgroundBiasEstimation.txt", stringsAsFactors=F)
  boxplot(tst$estimated - tst$simulated, range=0, ylab="estimated - simulated (degrees)")
  legend("topright", paste("N =", ntst), bty="n")

  # Result with ntst=1000, delta <2.5 degrees in 95% of the tests:
  # delta <- abs(tst$estimated - tst$simulated)
  # q <- quantile(delta, probs=c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm=T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
