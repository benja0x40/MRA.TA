# LIBRARIES ####################################################################

# -----------------------------------------------------------------------------.
# Load CRAN packages
# -----------------------------------------------------------------------------.
suppressPackageStartupMessages(library(methods))
suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(stringr))

# FUNCTIONS ####################################################################

# =============================================================================.
# Get command line arguments
# -----------------------------------------------------------------------------.
# Column 1: long option name
# Column 2: short option name 
# Column 3: 0=no argument, 1=required argument, 2=optional argument
# Column 4: data type (logical, integer, double, complex, character)
# Column 5: a brief description of the purpose of the option
# Column 6: default value
# -----------------------------------------------------------------------------.
processArgs <- function(
  script.args=NULL, auto.help=T, extra.help="", verbose=T
) {
  # ---------------------------------------------------------------------------.
  Args <- c()
  if(auto.help) {
    script.args <- rbind(
      c('help', 'h', 0, 'logical', "display usage informations", FALSE),
      script.args
    )
  }
  # ---------------------------------------------------------------------------.
  if(! is.null(script.args)) {
    getopt.specs <- rbind(script.args[,1:5])
    default.args <- script.args[,6]
    names(default.args) <- script.args[,1]
    arg.types <- script.args[,4] 
    names(arg.types) <- script.args[,1]
    Args <- getopt(spec = getopt.specs)
    Args$ARGS <- NULL
    for(pn in names(default.args)) {
      if(is.null(Args[[pn]])) {
        Args[[pn]] <- as(default.args[pn], Class=arg.types[pn])
      }
    }
    colnames(script.args) <- c(
      "name" ,"tag", "param", "type", "description", "default"
    )
  }
  # ---------------------------------------------------------------------------.
  if(auto.help) {
    if(Args$help) {
      cat("\n")
      cat(getopt(getopt.specs, usage=TRUE))
      cat("\n")
      cat(paste(paste("\t", extra.help, sep=""), collapse="\n"))
      cat("\n\n")
      q(save="no", status=1)
    }
  }
  # ---------------------------------------------------------------------------.
  n.args <- length(Args)
  if(verbose & n.args>0) {
    Args$script.name <- get_Rscript_filename()
    cat("\n", "[Rscript] ", Args$script.name, "\n")
    pn <- names(Args)[1:n.args]
    pn <- pn[pn!="help"]
    pv <- as.character(unlist(Args[pn]))
    pn <- str_pad(pn, max(nchar(pn))+2, side="right", pad=" ")
    cat(paste("\t", pn, "= ", pv, "\n", sep=""), "\n")
  }
  # ---------------------------------------------------------------------------.
  Args$script.args <- script.args
  Args
}
