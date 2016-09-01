

#' %% ~~ data name/kind ... ~~ File format definitions
#' 
#' %% ~~ A concise (1-5 lines) description of the dataset. ~~ A list of format
#' definitions for tab delimited text files. Available formats: \itemize{
#' \item\code{experimental.design} \item\code{mapped.probe.sequences}
#' \item\code{restriction.fragment} \item\code{restriction.site}
#' \item\code{nimblegen.ndf} \item\code{nimblegen.updated}
#' \item\code{nimblegen.pair} }
#' 
#' %% ~~ If necessary, more details than the __description__ above ~~
#' 
#' @name file.formats
#' @docType data
#' @format Each format definition provides the following informations:
#' \itemize{ \item\code{comment} : \code{character}, same as the
#' \code{comment.char} parameter of \code{read.table} function \item\code{skip}
#' : \code{integer}, same as the \code{skip} parameter of \code{read.table}
#' function \item\code{header} : \code{logical}, same as the \code{header}
#' parameter of \code{read.table} function \item\code{sep} : \code{character},
#' same as the \code{sep} parameter of \code{read.table} function
#' \item\code{table} : \code{data.frame} with the following fields \itemize{
#' \item\code{columns} : \code{character}, original column names when header is
#' present \item\code{classes} : \code{character}, R classes used to import
#' column values. Same as the \code{colClasses} parameter of \code{read.table}
#' function \item\code{extract} : \code{logical}, enable/disable the
#' importation of each column \item\code{name.as} : \code{character}, optional
#' renaming of extracted columns \item\code{description} : \code{character},
#' column description }
#' 
#' }
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' \code{\link{readData}}, \code{\link{verifyInputFiles}}
#' @references %% ~~ possibly secondary sources and usages ~~
#' @source %% ~~ reference to a publication or URL from which the data were
#' obtained ~~
#' @keywords datasets
#' @examples
#' 
#' data(file.formats)
#' 
NULL





#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_title(\"#1\")}",
#' "MRA.TA")\Sexpr{tools:::Rd_package_title("MRA.TA")}
#' 
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_description(\"#1\")}",
#' "MRA.TA")\Sexpr{tools:::Rd_package_description("MRA.TA")}
#' 
#' 
#' The DESCRIPTION file:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_DESCRIPTION(\"#1\")}",
#' "MRA.TA")\Sexpr{tools:::Rd_package_DESCRIPTION("MRA.TA")}
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_indices(\"#1\")}",
#' "MRA.TA")\Sexpr{tools:::Rd_package_indices("MRA.TA")} %% ~~ An overview of
#' how to use the package, including the most important functions ~~
#' 
#' @name MRA.TA-package
#' @aliases MRA.TA-package MRA.TA
#' @docType package
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "MRA.TA")\Sexpr{tools:::Rd_package_author("MRA.TA")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "MRA.TA")\Sexpr{tools:::Rd_package_maintainer("MRA.TA")}
#' @seealso %% ~~ Optional links to other man pages, e.g. ~~ %% ~~
#' \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
#' @references %% ~~ Literature or other references for background information
#' ~~
#' @keywords package
#' @examples
#' 
#' %% ~~ simple examples of the most important functions ~~
#' 
NULL



