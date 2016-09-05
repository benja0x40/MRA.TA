# =============================================================================.
#' kd-tree representation of (A,M) values
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
# -----------------------------------------------------------------------------.
#' @description
#' generate a kd-tree representation of (A,M) values for faster execution of the
#' background traversal algorithm (see normalizeArrayData).
# -----------------------------------------------------------------------------.
#' @param M
#' matrix of A and M values
#'
#' @param i
#' row indexes
#'
#' @param d
#'
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
kdtree <- function(M, i, d) {

  n <- length(i)
  node <- list()

  node$sg.tag <- FALSE
  node$sg.nbr <- n

  if(n==1) {
    node$sg.min <- M[i,]
    node$sg.max <- M[i,]
    node$sg.avg <- M[i,]
    node$sg.ids <- i
    return(node)
  }

  node$sg.min <- apply(M[i,], MARGIN=2, FUN=min)
  node$sg.max <- apply(M[i,], MARGIN=2, FUN=max)
  node$sg.avg <- apply(M[i,], MARGIN=2, FUN=mean)

  s.d <- abs(node$sg.max - node$sg.min)
  s.c <- as.integer(which.max(s.d))

  if( s.d[s.c] < d) {
    node$sg.ids <- i
  }
  else {
    s.m <- median(M[i, s.c])

    node$sg.tag <- TRUE
    node$sg.cut <- s.c
    if(s.m - min(M[i, s.c]) > max(M[i, s.c]) - s.m ) {
      node$sg.inf <- kdtree(M, i[ M[i, s.c] <  s.m ], d)
      node$sg.sup <- kdtree(M, i[ M[i, s.c] >= s.m ], d)
    }
    else {
      node$sg.inf <- kdtree(M, i[ M[i, s.c] <= s.m ], d)
      node$sg.sup <- kdtree(M, i[ M[i, s.c] >  s.m ], d)
    }
  }
  return(node)
}
# =============================================================================.
#' visualize the kd-tree structure
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
# -----------------------------------------------------------------------------.
#' @description
#' Utility function for visualizing the kd-tree structure of (A,M) values.
# -----------------------------------------------------------------------------.
#' @param node
# -----------------------------------------------------------------------------.
#' @return NULL
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
plot.kdtree.border <- function(node) {
  if(node$sg.tag) {
    plot.kdtree.border(node$sg.inf)
    plot.kdtree.border(node$sg.sup)
    rect(node$sg.min[1], node$sg.min[2], node$sg.max[1], node$sg.max[2], col=rgb(0,0,1,0.02), border=rgb(0,0,0,0.05))
  }
  else {
    rect(node$sg.min[1], node$sg.min[2], node$sg.max[1], node$sg.max[2], col=rgb(0,0,1,0.02), border=rgb(0,0,0,0.05))
  }
}
# =============================================================================.
#' Neighbor selection
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
# -----------------------------------------------------------------------------.
#' @description
#' Neighbor selection based on a pre-computed kd-tree structure of (A,M)
#' values.
# -----------------------------------------------------------------------------.
#' @param M
#'
#' @param v
#'
#' @param d
#'
#' @param node
#'
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
kdtree.getDN <- function(M, v, d, node) {

  n <- node$sg.nbr
  m <- length(v)

  if(node$sg.tag) {
    a <- c()
    b <- c()
    s.c <- node$sg.cut

    if((node$sg.inf$sg.max[s.c] + d) >= v[s.c]) {
      a <- kdtree.getDN(M, v, d, node$sg.inf)
    }
    if((node$sg.sup$sg.min[s.c] - d) <= v[s.c]) {
      b <- kdtree.getDN(M, v, d, node$sg.sup)
    }
    return(rbind(a, b))
  }

  i <- node$sg.ids
  A <- M[i,] - matrix(v, n, m, byrow=TRUE)
  A <- sqrt(apply(A*A, MARGIN=1, FUN=sum))
  match <- A <= d
  return(cbind(i[ match ], A[ match ]))
}
