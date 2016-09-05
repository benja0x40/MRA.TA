# =============================================================================.
#' Meanshift
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
#' @export
#' @seealso
#'    \link{},
# -----------------------------------------------------------------------------.
#' @description
#' Core function of the background traversal algorithm based on a pre-computed
#' kd-tree structure of (A,M) values. See normalizeArrayData function.
#' @details
# -----------------------------------------------------------------------------.
#' @param M
#' @param v
#' @param d
#' @param node
#' @param epsilon
#' @param col
#' @param nsteps
#' @param plots
# -----------------------------------------------------------------------------.
#' @return
# -----------------------------------------------------------------------------.
#' @examples
# -----------------------------------------------------------------------------.
mean.shift <- function(M, v, d, node, epsilon=0.001, col=rgb(1,0.5,0,1), nsteps=0, plots=FALSE) {

  i <- rep(FALSE, nrow(M))
  m <- length(v)
  moves <- v
  alpha <- c()

  while(TRUE) {

    IW <- kdtree.getDN(M, v, d, node)

    n <- nrow(IW)

    if(length(n)==0)
      return(NULL)

    i[IW[,1]] <- TRUE

    A <- M[IW[,1],] - matrix(v, n, m, byrow=TRUE)
    A <- apply(A, MARGIN=2, FUN=mean)

    if(plots) segments(v[1], v[2], (v+A)[1], (v+A)[2], col=col,lwd=2)

    v <- v + A
    moves <- rbind(moves, v)

    a <- asin(A[2]/sqrt(A[1]*A[1]+A[2]*A[2]))
    alpha <- c(alpha, a)

    if(length(alpha)>nsteps & sqrt(sum(A*A)) < epsilon) break
  }

  list(i=i, centroid=v, alpha=alpha, moves=moves)
}
