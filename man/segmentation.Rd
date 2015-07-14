\name{segmentation}
\alias{segmentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Multi-resolution domain segmentation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
segmentation(Yi, name, wmax = 0, wmin = 3, gamma = 0.01, Tw = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Yi}{
    %%     ~~Describe \code{Yi} here~~
    vector of statistical scores produced by the \code{enrichmentScore}
    function.
  }
  \item{name}{
    %%     ~~Describe \code{name} here~~
    prefix used to generate output file names.
  }
  \item{wmax}{
    %%     ~~Describe \code{wmax} here~~
    maximal window size to be scanned in number of Yi values.
    The default is \code{length(Yi)}.
  }
  \item{wmin}{
    %%     ~~Describe \code{wmin} here~~
    minimum size of segmented domains in number of Yi values.
    The default is 3 consecutive Yi values.
  }
  \item{gamma}{
    %%     ~~Describe \code{gamma} here~~
    stringency factor ranging from ]0;1]. Default value is 1.0 and
    ower values correspond to increased stringency.
  }
  \item{Tw}{
    %%     ~~Describe \code{Tw} here~~
    a log value ranging from ]\code{-Inf};0] and indicating the minimal statistical
    scores (e.g. pseudo p-values from combined rank-based scores) accepted for
    segmented domains.
    Default value is 0 and equivalent to an absence of cutoff while lower values
    correspond to more stringent cutoffs.
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
  List object with the following attributes:
  \item{time.optimization}{computation time for the optimization procedure}
  \item{time.MRT.Analysis}{computation time for the domain fusion procedure}
  \item{n.segments}{total number of locally optimal segments after optimization}
  \item{n.domains}{total number of domains after domain fusion}
  \item{n.maxresolution}{number of maximum resolution domains}
  \item{n.maxscale}{number of maximum scale domains}
  \item{file.segments}{result file for the optimization procedure}
  \item{file.domains}{result file for the domain fusion procedure}
  \item{file.maxresolution}{result file for maximum resolution domains}
  \item{file.maxscale}{result file for maximum scale domains}
  \item{Pw.min}{
    best statistical score (pseudo p-values from combined rank-based scores)
    for each scanned window size
  }
}
\section{Output files}{
  \itemize{
    \item{
      Common structure
      
      All output files are tab delimited text files including two header lines.
      \preformatted{
        line 1          = analysis parameters
        line 2          = column names
        following lines = data table
      }
    }
    \item{
    \code{file.segments}
    
    This file contains the optimal segment data,
    resulting from the the optimization procedure,
    and defined by columns (i, w, Piw) as follows:
      \preformatted{
        i   = index of the last Yi value within the optimal segment
        w   = size in number of consecutive Yi values
        Piw = statistical score
      }
    }
    \item{
      \code{file.domains}, \code{file.maxresolution}, \code{file.maxscale}
      
      These files contain the segmented domain data,
      resulting from the domain fusion procedure,
      and defined by columns (id, container, start, end, wmin, wmax, P, i, w)
      as follows:
      \preformatted{
        id        = unique identifier for each domain
        container = identifier of parent domain, 0 meaning no parent
        start     = index of the first Yi value within domain
        end       = index of the last Yi value within domain
        wmin      = minimum size of included optimal segments, in number of Yi values
        wmax      = maximum size of included optimal segments, in number of Yi values
        P         = statistical score of the locally optimal segment
        i         = index of the last Yi value within the locally optimal segment
        w         = size of the locally optimal segment in number of consecutive Yi values
      }
    }
  }
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
  \code{\link{calc.Qi}},
  \code{\link{enrichmentScore}},
  \code{\link{domainogram}},
  \code{\link{plotOptimalSegments}},
  \code{\link{plotDomains}}
}
\examples{
  # Simulate enrichment signal
  n <- 2000
  Mi <- rep(0, n)
  Mi <- Mi + dnorm(1:n, 2.5*n/20,  n/40) + dnorm(1:n, 4*n/20,  50)
  Mi <- Mi + 4 * dnorm(1:n, 5*n/10, n/10)
  Mi <- Mi + dnorm(1:n, 16*n/20,  n/40) + dnorm(1:n, 17.5*n/20,  50)
  Mi <- (Mi/max(Mi))^4 + rnorm(n)/4
  
  # Compute enrichment scores
  Yi <- enrichmentScore(Mi)
  
  # Multi-resolution segmentation
  seg.c <- segmentation(Yi, name="MRA_demo", wmin=20)
  
  # Load segmentation results
  opts <- read.delim(seg.c$file.segments, stringsAsFactors=F, skip=1)
  doms <- read.delim(seg.c$file.domains, stringsAsFactors=F)
  doms.mr <- read.delim(seg.c$file.maxresolution, stringsAsFactors=F)
  doms.ms <- read.delim(seg.c$file.maxscale, stringsAsFactors=F)
  
  # Visualization coordinates
  x.s <- 1:n - 0.5
  x.e <- 1:n + 0.5
  w2y <- wSize2yAxis(n, logscale=T)
  
  layout(matrix(1:2, 2, 1), heights=c(3,1)/4)
  par(mar=c(3, 4, 1, 2)) # default bottom, left, top, right = c(5, 4, 4, 2)
  
  # Plot domainogram
  domainogram(Yi, x.s, x.e, w2y)
  plot(Mi, type='l')
  
  # Visualize segmentation results
  plotOptimalSegments(opts, x.s, x.e, w2y, col="black")
  plot(Mi, type='l')
  
  # Visualize multi-resolution domains
  plotDomains(doms, x.s, x.e, w2y, col=rgb(0,0,0,0.5), border=rgb(0,0,0,0))
  # Visualize max. resolution and max. scale domains
  plotDomains(doms.mr, x.s, x.e, w2y, col=rgb(0,1,0,0.5), border=rgb(0,0,0,0), lwd=2, lty=1, add=T)
  plotDomains(doms.ms, x.s, x.e, w2y, col=rgb(1,0,0,0.5), border=rgb(0,0,0,0), lwd=2, lty=1, add=T)
  legend("topright", c("Max. resolution", "Max. scale", "Both"), fill=c("green", "red", "chocolate"), bty='n')
  plot(Mi, type='l')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
