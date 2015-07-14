\name{computeRestrictionMap}
\alias{computeRestrictionMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute a restriction map 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeRestrictionMap(genome.seq, enzyme.motif, output.file, seqlist = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genome.seq}{
    %%     ~~Describe \code{genome.seq} here~~
    \code{DNAStringSet} object containing the genome sequence
  }
  \item{enzyme.motif}{
    %%     ~~Describe \code{enzyme.motif} here~~
    \code{character} value defining the enzyme cut sites as a DNA sequence
  }
  \item{output.file}{
    %%     ~~Describe \code{name} here~~
    \code{character} prefix used to generate output file names.
  }
  \item{seqlist}{
    %%     ~~Describe \code{seqlist} here~~
    indices of chromosomes to be considered, default = all
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
  \item{sites}{file name with the restriction sites data}
  \item{fragments}{file name with restriction fragments data}
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
  \code{\link{matchFragmentsToProbes}},
  \code{\link{matchProbesToFragments}},
  \code{\link{DNAStringSet}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
