matrix_to_list_of_rows = function(x)
  lapply(seq_len(nrow(x)), function(i) x[i,])

symmetrize_matrix = function(X) (X + t(X))/2

#' @export
printf = function(...) cat(sprintf(...))

#' @export
tic = function() proc.time()["elapsed"]

#' @export
toc = function(t0) proc.time()["elapsed"] - t0

#' @export
check_01_cols = function(X) {
  any(Matrix::colSums(X != 0) <= 1)
}


#' @export
remove_sparse_rowcol = function(X, thersh=1)  {
  prev_dim = c(0,0)
  new_dim = dim(X)
  while(sum(abs(new_dim - prev_dim)) > 0) {
    cat('.')
    X <- X[Matrix::rowSums(X > 0) > thersh,]
    X <- X[, Matrix::colSums(X > 0) > thersh]
    prev_dim = new_dim
    new_dim = prod(dim(X))
  }
  cat('\n')
  return(X)
}

#' @export
remove_sparse_cols = function(X, thersh=1)  {
    X[ , Matrix::colSums(X > 0) > thersh]
}

#' @export
imagesc <- function(X, remove_diag=T, scale=T) {
  X <- Matrix::Matrix(X)
  if (remove_diag) diag(X) = 0
  if (scale)
    Matrix::image(cbind(c(c(1,-1), rep(0,nrow(X)-2)), X),
                  xlab="",ylab="", sub=NA, border.col=NA, xlim=c(2,ncol(X)+1))
  else
    Matrix::image(X,xlab="",ylab="", sub=NA, border.col=NA)
}

#' @export
order_cols_by_freq = function(X, trunc = ncol(X), decreasing = T) {
  idx <- order(Matrix::colMeans(X), decreasing = decreasing)
  d <- min(trunc, nrow(X))
  X = X[ , idx[1:d] ]
  return(X)
}


