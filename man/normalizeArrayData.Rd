\name{normalizeArrayData}
\alias{normalizeArrayData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Normalize tiling array data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normalizeArrayData(A, M, smoothness = 0.08, epsilon = 0.01, nsteps = 11, name = "Test", plots = TRUE, lowess = T, lowess.f = 0.2, lowess.mad = 0, lowess.iter = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
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
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{plots}{
%%     ~~Describe \code{plots} here~~
}
  \item{lowess}{
%%     ~~Describe \code{lowess} here~~
}
  \item{lowess.f}{
%%     ~~Describe \code{lowess.f} here~~
}
  \item{lowess.mad}{
%%     ~~Describe \code{lowess.mad} here~~
}
  \item{lowess.iter}{
%%     ~~Describe \code{lowess.iter} here~~
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
  \code{\link{backgroundBiasCorrection}},
  \code{\link{lowessCorrection}}
}
\examples{
  # Load array data and apply background bias estimation+correction, followed by
  # lowess normalization
  
  data(WT.4C.Fab7.dm6)

  # 1. Combined procedure ------------------------------------------------------
  
  # Raw A and M values
  A <- (log2(r1.4C$PM) + log2(r1.ct$PM))/2
  M <- (log2(r1.4C$PM) - log2(r1.ct$PM))
  
  # Normalized A and M values
  res <- normalizeArrayData(
    A, M, name="4C_norm", plots=TRUE
  )
  A <- res$A; M <- res$M
  
  # 2. Equivalent step by step procedure ---------------------------------------
  
  # Raw A and M values
  A <- (log2(r1.4C$PM) + log2(r1.ct$PM))/2
  M <- (log2(r1.4C$PM) - log2(r1.ct$PM))

  # Estimate background bias
  bb.r1 <- backgroundBiasEstimation(A, M, plots = T)
  
  # Correct background bias
  res <- backgroundBiasCorrection(A, M, theta=bb.r1)
  A <- res$x; M <- res$y
  
  # Apply lowess normalization
  res <- lowessCorrection(A, M, lowess.f=0.2, plots = T)
  A <- res$x; M <- res$y
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
