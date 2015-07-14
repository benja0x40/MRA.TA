\name{backgroundBiasCorrection}
\alias{backgroundBiasCorrection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Background bias correction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Correct for a background bias as estimated by the \code{backgroundBiasEstimation} function.
}
\usage{
backgroundBiasCorrection(x, y, theta, AM.scale.compensation = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{AM.scale.compensation}{
%%     ~~Describe \code{AM.scale.compensation} here~~
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
  \code{\link{backgroundBiasEstimation}},
  \code{\link{normalizeArrayData}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
