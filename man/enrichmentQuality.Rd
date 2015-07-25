\name{enrichmentQuality}
\alias{enrichmentQuality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Best and worst enrichment levels
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Identify probes most likely associated with best and worst enrichment levels.
}
\usage{
enrichmentQuality(x, y, q.best = 0.015, q.worst = 0.005, plots = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
  }
  \item{y}{
    %%     ~~Describe \code{y} here~~
  }
  \item{q.best}{
    %%     ~~Describe \code{q.best} here~~
  }
  \item{q.worst}{
    %%     ~~Describe \code{q.worst} here~~
  }
  \item{plots}{
    %%     ~~Describe \code{plots} here~~
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
\code{\link{plotProbeDistanceControls}}
}
\examples{
  x <- 2^pmin(rlnorm(100000, 1, 0.5), 10)
  y <- 2^pmin(rlnorm(100000, 0.5, 1), 10)

  chk <- enrichmentQuality(x, y, plots=T)

  # Visualize thresholds on enrichment quality scores
  clr <- with(
    chk, rgb(is.worst, is.best, 0, ifelse(is.worst | is.best, 0.5, 0.1))
  )
  with(chk, plot(w.score, b.score, pch='.', col=clr))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
