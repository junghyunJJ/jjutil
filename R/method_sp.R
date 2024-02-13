convert_dataframe <- function(dat, preview = TRUE) {
  dat <- as.data.frame(dat)
  rownames(dat) <- dat[, 1]
  dat <- dat[, -1]
  if (preview) {
    print(dat[1:5, 1:5])
  }
  return(dat)
}

linear_kernel <- function(X) {
  K <- tcrossprod(as.matrix(X)) # K = np.dot(X, X.T)
  return(K / max(K))            # return K / K.max()
}
