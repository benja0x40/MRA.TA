\name{wSize2yAxis}
\alias{wSize2yAxis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Continuous coordinates for domainogram visualization (y axis)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compute coordinates for the domainogram visualization of a range of window sizes
(y axis).
}
\usage{
wSize2yAxis(wmax, logscale = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wmax}{
%%     ~~Describe \code{wmax} here~~
integer indicating the maximum window size to be represented
}
  \item{logscale}{
%%     ~~Describe \code{logscale} here~~
logical indicating if window sizes should be represented according to a log2
scale instead of a linear scale (default)
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
  \item{logscale}{value of the \code{logscale} parameter}
  \item{ylab}{text to be used for the y axis legend}
  \item{y.t}{top coordinate for each window size}
  \item{y.b}{bottom coordinate for each window size}
  \item{ylim}{range of the visualization coordinates (y axis)}
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
\code{\link{domainogram}},
\code{\link{visualizationCoordinates}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
