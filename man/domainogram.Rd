\name{domainogram}
\alias{domainogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Domainogram visualization
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Plot a domainogram from a sequence of log-transformed rank-based scores.

  This function allows to use a linear or log2 scale for window sizes
  (y axis).
}
\usage{
domainogram(Yi, x.s, x.e, w2y, xlim = NULL, wlim = NULL, colfun = makeColors, xlab = "", add = F, logskip = 1.01, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Yi}{
    %%     ~~Describe \code{Yi} here~~
    vector of statistical scores produced by the \code{enrichmentScore}
    function
  }
  \item{x.s}{
    %%     ~~Describe \code{x.s} here~~
    vector of start coordinates (x axis) defining the represented genomic intervals
  }
  \item{x.e}{
    %%     ~~Describe \code{x.e} here~~
    vector of end coordinates (x axis) defining the represented genomic intervals
  }
  \item{w2y}{
    %%     ~~Describe \code{w2y} here~~
    coordinate mapping for window sizes (y axis) produced by the
    \code{wSize2yAxis} function
  }
  \item{xlim}{
    %%     ~~Describe \code{xlim} here~~
    range of the represented genomic region = c(xmin, xmax)
  }
  \item{wlim}{
    %%     ~~Describe \code{wlim} here~~
    range of represented window sizes = c(wmin, wmax)
  }
  \item{colfun}{
    %%     ~~Describe \code{colfun} here~~
    coloring function, uses \code{makeColors} as the default coloring method.
  }
  \item{xlab}{
    %%     ~~Describe \code{xlab} here~~
    legend of the x axis
  }
  \item{add}{
    %%     ~~Describe \code{add} here~~
    \code{logical}, if set to TRUE indicates to overlay the domainogram on an
    existing and consistent multi-resolution plot
  }
  \item{logskip}{
    %%     ~~Describe \code{logskip} here~~
    \code{numeric} coefficient for skipping the representation of window sizes
    when using the log2 scale on the y axis. A value of 1.0 forces to
    represent all window size. A value of 2 would double the window sizes
    represented at each iteration of the rendering.
  }
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
    optional parameters forwarded to the \code{colfun} function. When using
    the default \code{makeColors} function, such parameters are
    \code{thresholds}, \code{colors}, \code{background} and \code{overflow}
    (see documentation for details).
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  The domainogram is a multi-resolution representation of scan statistics 
  introduced by de Wit et al., 2008 (see references below). It is based on
  R.A. Fisher's method for the combination of p-values applied to 
  non-parametric scores.
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
de Wit E., Braunschweig U., Greil F., Bussemaker H.J., and van Steensel B.,
2008. Global chromatin domain organization of the Drosophila genome.
PLoS genetics 4: e1000045. \url{http://www.ncbi.nlm.nih.gov/pubmed/18369463}
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
  \code{\link{calc.Qi}},
  \code{\link{enrichmentScore}},
  \code{\link{visualizationCoordinates}},
  \code{\link{wSize2yAxis}},
  \code{\link{makeColors}},
  \code{\link{segmentation}},
  \code{\link{plotOptimalSegments}},
  \code{\link{plotDomains}}
}
\examples{
  # Simulate enrichment signal
  n <- 500
  Mi <- cos(2*pi*2.5*(-n:n)/n) + rnorm(2*n+1)
  n <- length(Mi)
  
  # Compute enrichment scores
  Yi <- enrichmentScore(Mi)
  
  # Visualization coordinates
  x.s <- 1:n - 0.5
  x.e <- 1:n + 0.5
  w2y <- wSize2yAxis(n)
  
  layout(matrix(1:2, 2, 1), heights=c(3,1)/4)
  par(mar=c(3, 4, 1, 2)) # default bottom, left, top, right = c(5, 4, 4, 2)
  
  # Plot domainogram
  domainogram(Yi, x.s, x.e, w2y)
  plot(Mi, type='l')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
