\name{file.formats}
\alias{file.formats}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
File format definitions
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
A list of format definitions for tab delimited text files.
Available formats:
  \itemize{
    \item{\code{experimental.design}}{ }
    \item{\code{mapped.probe.sequences}}{ }
    \item{\code{restriction.fragment}}{ }
    \item{\code{restriction.site}}{ }
    \item{\code{nimblegen.ndf}}{ }
    \item{\code{nimblegen.updated}}{ }
    \item{\code{nimblegen.pair}}{ }
  }
}
\usage{data("file.formats")}
\format{
Each format definition provides the following informations:
  \itemize{
    \item{\code{comment}}{ : \code{character}, same as the \code{comment.char} parameter of \code{read.table} function }
    \item{\code{skip}}{ : \code{integer}, same as the \code{skip} parameter of \code{read.table} function }
    \item{\code{header}}{ : \code{logical}, same as the \code{header} parameter of \code{read.table} function }
    \item{\code{sep}}{ : \code{character}, same as the \code{sep} parameter of \code{read.table} function }
    \item{\code{table}}{ : \code{data.frame} with the following fields
      \itemize{
        \item{\code{columns}}{ : \code{character}, original column names when header is present }
        \item{\code{classes}}{ : \code{character}, R classes used to import column values. Same as the \code{colClasses} parameter of \code{read.table} function }
        \item{\code{extract}}{ : \code{logical}, enable/disable the importation of each column }
        \item{\code{name.as}}{ : \code{character}, optional renaming of extracted columns }
        \item{\code{description}}{ : \code{character}, column description }
      }
    }
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{readData}},
\code{\link{verifyInputFiles}}
}
\examples{
data(file.formats)
}
\keyword{datasets}
