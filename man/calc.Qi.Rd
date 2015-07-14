\name{calc.Qi}
\alias{calc.Qi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Rank based scores
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Transform numeric values into rank-based scores.
}
\usage{
calc.Qi(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
vector of numeric values
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
  calc.Qi(x) = (rank(x) - 1) / n\cr
  where n is the length of x and ranking is performed in decreasing order:
  \itemize{
    \item rank(min(x)) = n
    \item rank(max(x)) = 1
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
\code{\link{enrichmentScore}},
\code{\link{domainogram}},
\code{\link{segmentation}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
