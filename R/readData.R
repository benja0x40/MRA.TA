# =============================================================================.
#' Read tabulated text files
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export readData
#' @seealso
#'    \link{verifyInputFiles}
# -----------------------------------------------------------------------------.
#' @description
#' Read tab delimited text files with specific format definition
#'
#' @details
#' List of available format definitions for tab delimited text files:
#' \itemize{
#'   \item \code{experimental.design}
#'   \item \code{mapped.probe.sequences}
#'   \item \code{restriction.fragment}
#'   \item \code{restriction.site}
#'   \item \code{nimblegen.ndf}
#'   \item \code{nimblegen.updated}
#'   \item \code{nimblegen.pair}
#' }
#'
#' Each format definition provides the following informations:
#' \itemize{
#'   \item \code{comment} : \code{character}, same as the \code{comment.char} parameter of \code{read.table} function
#'   \item\code{skip}     : \code{integer}, same as the \code{skip} parameter of \code{read.table} function
#'   \item\code{header}   : \code{logical}, same as the \code{header} parameter of \code{read.table} function
#'   \item\code{sep}      : \code{character}, same as the \code{sep} parameter of \code{read.table} function
#'   \item\code{table}    : \code{data.frame} with the following fields
#'   \itemize{
#'     \item\code{columns}  : \code{character}, original column names when header is present
#'     \item\code{classes}  : \code{character}, R classes used to import column values. Same as the \code{colClasses} parameter of \code{read.table} function
#'     \item\code{extract}  : \code{logical}, enable/disable the importation of each column
#'     \item\code{name.as}  : \code{character}, optional renaming of extracted columns
#'     \item\code{description} : \code{character}, column description
#'   }
#' }
# -----------------------------------------------------------------------------.
#' @param fpath
#'
#' @param fform
#' file format definition (see details below).
#'
# -----------------------------------------------------------------------------.
#' @return
#' readData returns a data.frame
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
readData <- function(fpath, fform) {
  colClasses <- c("NULL", "integer", "numeric", "character", "logical")
  if(is.null(fform$comment)) fform$comment <- "#"
  if(is.null(fform$skip))    fform$skip    <- 0
  if(is.null(fform$header))  fform$header  <- T
  if(is.null(fform$sep))     fform$sep     <- "\t"
  if(is.null(fform$table)) {
    stop("Missing table format", call. = F)
  }
  if(! "columns" %in% colnames(fform$table)) {
    # *** WARNING/TODO *** #
  }
  if(! "classes" %in% colnames(fform$table)) {
    # *** WARNING/TODO *** #
  }
  if(sum(! fform$table$classes %in% colClasses)>0) {
    stop("Incorrect column classes", call. = F)
  }
  if(is.null(fform$table$extract)) {
    fform$table$extract <- rep(T, nrow(fform$table))
  }
  if(is.null(fform$table$rename.as)) {
    fform$table$rename.as <- rep(NA, nrow(fform$table))
  }
  xtr <- fform$table$extract & (fform$table$classes != "NULL")
  cls <- fform$table$classes
  cls[! xtr] <- "NULL"
  lbl <- fform$table$columns
  idx <- which(! is.na(fform$table$rename.as))
  if(length(idx)>0) lbl[idx] <- fform$table$rename.as[idx]
  # Load table data
  tbl <- read.table(
    fpath,
    comment.char=fform$comment,  skip=fform$skip, header=fform$header,
    sep=fform$sep, colClasses=cls
  )
  if(! is.null(lbl)) {
    colnames(tbl) <- lbl[xtr]
  }
  tbl
}
