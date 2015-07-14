\name{visualizationCoordinates}
\alias{visualizationCoordinates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Continuous coordinates for domainogram visualization (x axis)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compute coordinates for the domainogram visualization of a range of genomic
positions (x axis).
}
\usage{
visualizationCoordinates(start, end)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
%%     ~~Describe \code{start} here~~
vector defining the start positions of possibily disjoint genomic
intervals
}
  \item{end}{
%%     ~~Describe \code{end} here~~
vector defining the end positions of possibily disjoint genomic
intervals
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
\item{start }{start positions for domainogram visualization}
\item{end }{end positions for domainogram visualization}
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
\code{\link{wSize2yAxis}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
