# LIBRARIES ####################################################################

# -----------------------------------------------------------------------------.
# Load CRAN packages
# -----------------------------------------------------------------------------.
suppressPackageStartupMessages(library("stringr"))

# FUNCTIONS ####################################################################

# =============================================================================.
# Verify input files
# -----------------------------------------------------------------------------.
verifyInputFiles <- function(path.list, columns=NULL, sep="\t", comments="#") {
  idx <- which(! file.exists(path.list))
  if(length(idx)>0) {
    msg <- paste(unique(path.list[idx]), collapse="\n\t")
    msg <- paste("File(s) not found:\n\t", msg, "\n", sep="")
    stop(msg, call.=F)
  }
  if(! is.null(columns)) {
    rex <- paste("^", comments, sep="")
    for(fpath in path.list) {
      nbr <- 0
      tst <- comments
      while(grepl(rex, tst)) {
        nbr <- nbr+1
        tst <- readLines(fpath, n=nbr)[nbr]
      }
      tst <- unlist(str_split(tst, sep))
      idx <- which(! columns %in% tst)
      if(length(idx)>0) {
        msg <- paste(columns[idx], collapse="\n\t")
        msg <- paste("Missing columns in ", fpath, "\n\t", msg, "\n", sep="")
        stop(msg, call.=F)
      }
    }
  }
}
# =============================================================================.
# 
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

# FORMAT DEFINITIONS ###########################################################

file.formats <- list()

# Experimental Designs =========================================================

# -----------------------------------------------------------------------------.
# Experimental Design
file.formats$experimental.design <- list(
  comment ="#",
  skip    = 0,
  header  = T,
  sep     = "\t",
  table = data.frame(
    columns = c(
      "ID",          # character 
      "REPLICATE",   # character
      "LIBRARY",     # character
      "SAMPLE_PATH", # character
      "DESIGN_PATH"  # character
    ),
    classes = c(
      "character",    # ID
      "character",    # REPLICATE
      "character",    # LIBRARY
      "character",    # SAMPLE_PATH
      "character"     # DESIGN_PATH
    ),
    extract = rep(T, 5),
    name.as = rep(NA, 5),
    description = rep("", 5),
    stringsAsFactors=F
  )
)

# Mapped Sequences =============================================================

# -----------------------------------------------------------------------------.
# Bowtie alignment from the fasta sequences of array probes
file.formats$mapped.probe.sequences <- list(
  comment ="",
  skip    = 0,
  header  = F,
  sep     = "\t",
  table = data.frame(
    classes = c(
      "integer",    # 1: read index
      "character",  # 2: PROBE_ID
      "character",  # 3: strand
      "character",  # 4: chromosome/seq_id
      "integer",    # 5: position
      "character",  # 6: PROBE_SEQUENCE
      "character",  # 7: matching nucleotides
      "integer",    # 8: ???
      "character"   # 9: ???
    ),
    extract = rep(T, 9),
    name.as = rep(NA, 9),
    description = rep("", 9),
    stringsAsFactors=F
  )
)

# Restriction Maps =============================================================

# -----------------------------------------------------------------------------.
# Restriction Fragments
file.formats$restriction.fragment <- list(
  comment ="#",
  skip    = 0,
  header  = T,
  sep     = "\t",
  table = data.frame(
    columns = c(
      "RFID",      # character 
      "CHR",       # character
      "RS5.START", # integer
      "RS5.END",   # integer
      "RS3.START", # integer
      "RS3.END"    # integer
    ),
    classes = c(
      "character", # RFID
      "character", # CHR
      "numeric",   # RS5.START
      "numeric",   # RS5.END
      "numeric",   # RS3.START
      "numeric"    # RS3.END
    ),
    extract = rep(T, 6),
    name.as = rep(NA, 6),
    description = rep("", 6),
    stringsAsFactors=F
  )
)
# -----------------------------------------------------------------------------.
# Restriction Sites
file.formats$restriction.site <- list(
  comment ="#",
  skip    = 0,
  header  = T,
  sep     = "\t",
  table = data.frame(
    columns = c(
      "CHR",   # character
      "START", # numeric
      "END"    # numeric
    ),
    classes = c(
      "character", # CHR
      "numeric",   # START
      "numeric"    # END
    ),
    extract = rep(T, 3),
    name.as = rep(NA, 3),
    description = rep("", 3),
    stringsAsFactors=F
  )
)

# Array Design =================================================================

# -----------------------------------------------------------------------------.
# Nimblegen .ndf
file.formats$nimblegen.ndf <- list(
  comment ="#",
  skip    = 0,
  header  = T,
  sep     = "\t",
  table = data.frame(
    columns = c(
      "PROBE_DESIGN_ID",    # character 
      "CONTAINER",          # character
      "DESIGN_NOTE",        # character
      "SELECTION_CRITERIA", # character
      "SEQ_ID",             # character
      "PROBE_SEQUENCE",     # character
      "MISMATCH",           # integer
      "MATCH_INDEX",        # integer
      "FEATURE_ID",         # integer
      "ROW_NUM",            # integer
      "COL_NUM",            # integer
      "PROBE_CLASS",        # character
      "PROBE_ID",           # character
      "POSITION",           # integer
      "DESIGN_ID",          # integer
      "X",                  # integer
      "Y"                   # integer
    ),
    classes = c(
      "character",    # PROBE_DESIGN_ID
      "character",    # CONTAINER
      "character",    # DESIGN_NOTE
      "character",    # SELECTION_CRITERIA
      "character",    # SEQ_ID
      "character",    # PROBE_SEQUENCE
      "integer",      # MISMATCH
      "integer",      # MATCH_INDEX
      "integer",      # FEATURE_ID
      "integer",      # ROW_NUM
      "integer",      # COL_NUM
      "character",    # PROBE_CLASS
      "character",    # PROBE_ID
      "integer",      # POSITION
      "integer",      # DESIGN_ID
      "integer",      # X
      "integer"       # Y
    ),
    extract = rep(T, 17),
    name.as = rep(NA, 17),
    description = rep("", 17),
    stringsAsFactors=F
  )
)

# -----------------------------------------------------------------------------.
# Nimblegen array design updated (remapped to genome assembly)
file.formats$nimblegen.updated <- list(
  comment ="#",
  skip    = 0,
  header  = T,
  sep     = "\t",
  table = data.frame(
    columns = c(
      "PROBE_ID",    # character
      "PROBE_CLASS", # character
      "CONTAINER",   # character
      "X",           # integer
      "Y",           # integer
      "SEQ_ID",      # character
      "CHR",         # character
      "STRAND",      # character
      "START",       # integer
      "END",         # integer
      "LENGTH"       # integer   = length of the probe
    ),
    classes = c(
      "character", # PROBE_ID
      "character", # PROBE_CLASS
      "character", # CONTAINER
      "integer",   # X
      "integer",   # Y
      "character", # SEQ_ID
      "character", # CHR
      "character", # STRAND
      "integer",   # START
      "integer",   # END
      "integer"    # LENGTH
    ),
    extract = rep(T, 11),
    name.as = rep(NA, 11),
    description = rep("", 11),
    stringsAsFactors=F
  )
)

# Array Data ===================================================================

# -----------------------------------------------------------------------------.
# Nimblegen hybridization scans (pair files)
file.formats$nimblegen.pair <- list(
  comment ="#",
  skip    = 0,
  header  = T,
  sep     = "\t",
  table = data.frame(
    columns = c(
      "IMAGE_ID",         # character
      "GENE_EXPR_OPTION", # character
      "SEQ_ID",           # character
      "PROBE_ID",         # character
      "POSITION",         # numeric
      "X",                # integer
      "Y",                # integer
      "MATCH_INDEX",      # integer
      "SEQ_URL",          # character
      "PM",               # numeric
      "MM"                # numeric
    ),
    classes = c(
      "character", # IMAGE_ID
      "character", # GENE_EXPR_OPTION
      "character", # SEQ_ID
      "character", # PROBE_ID
      "numeric",   # POSITION
      "integer",   # X
      "integer",   # Y
      "integer",   # MATCH_INDEX
      "character", # SEQ_URL
      "numeric",   # PM
      "numeric"    # MM
    ),
    extract = c(F, F, F, T, F, T, T, F, F, T, F),
    name.as = rep(NA, 11),
    description = rep("", 11),
    stringsAsFactors=F
  )
)

# # ------------------------------------------------------------------------------
# # Read micro array data from the Agilent platform
# readAgilent <- function(filename) {
#   col.classes <- c(
#     rep("NULL",2),
#     rep("integer",2),
#     rep("NULL",4),
#     "character",
#     rep("NULL",2),
#     "character",
#     "NULL",
#     "character",
#     rep("NULL",3),
#     rep("double",3),
#     rep("NULL",4),
#     rep("double",4),
#     rep("NULL",85)
#   )
#   target <- col.classes != "NULL"
#   names <- readLines(filename, 10)
#   names <- names[10]
#   names <- unlist(strsplit(names,split="\t"))
#   names <- names[target]
#   data <- cgh.data1 <- read.table(
#     file=filename,
#     skip=10,
#     header=FALSE,
#     sep="\t",
#     quote="",
#     dec=".",
#     colClasses=col.classes
#   )
#   colnames(data) <- names
#   chr   <- rep("",nrow(data))
#   start <- rep(0,nrow(data))
#   end   <- rep(0,nrow(data))
#   pb <- txtProgressBar(min = 1, max = nrow(data), style = 3, char="|")
#   for(i in 1:nrow(data)) {
#     tmp <- unlist(strsplit(data[i,"SystematicName"],":"))
#     chr[i] <- tmp[1]
#     tmp <- unlist(strsplit(tmp[2],"-"))
#     start[i] <- as.integer(tmp[1])
#     end[i] <- as.integer(tmp[2])
#     setTxtProgressBar(pb, i)
#   }
#   close(pb)
#   data$original.chr <- chr
#   data$original.start <- start
#   data$original.end <- end
#   
#   data
# }
# # -----------------------------------------------------------------------------.
# agl.file <- "Raw_Data/GSE23166_RAW/GSM570307_WT_Abd-B_rep1a.txt.gz"
# agl.head <- readLines(agl.file, 12)
# agl.type <- unlist(strsplit(agl.head[9], split="\t"))
# agl.cols <- unlist(strsplit(agl.head[10], split="\t"))
# agl.ncol <- length(agl.cols)
# # agl.type <- agl.type[2:agl.ncol]
# # agl.cols <- agl.cols[2:agl.ncol]
# agl.type <- gsub("TYPE|text", "character", agl.type)
# agl.type <- gsub("float|double", "numeric", agl.type)
# agl.type <- gsub("boolean", "integer", agl.type)
# file.formats$agilent <- list(
#   table = data.frame(
#     columns = agl.cols,
#     classes = agl.type,
#     extract = c(rep(F, 2), rep(T, 2), rep(F, 2), T, F, rep(T, agl.ncol - 8)),
#     name.as = c(rep(NA, 2), c("X", "Y"), rep(NA, agl.ncol - 4)),
#     stringsAsFactors=F
#   ),
#   skip = 10,
#   header=F,
#   comment="#"
# )
# 
# # Agilent hybridization scans (pair files)
# file.formats$agilent <- data.frame(
#   colnames = c(
#   ),
#   colClasses = c(
#     rep("NULL",2),
#     rep("integer",2),
#     rep("NULL",4),
#     "character",
#     rep("NULL",2),
#     "character",
#     "NULL",
#     "character",
#     rep("NULL",3),
#     rep("double",3),
#     rep("NULL",4),
#     rep("double",4),
#     rep("NULL",85)
#   ),
#   stringsAsFactors=F
# )
# write.table(file.formats$experimental.design$table, "R/Table_Formats/experimental.design", sep="\t", col.names=T, row.names=F, quote=F)
