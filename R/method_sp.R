convert_dataframe <- function(dat, preview = FALSE) {
  dat <- as.data.frame(dat)
  rownames(dat) <- dat[, 1]
  dat <- dat[, -1, drop = FALSE]

  if (preview) {
    if(ncol(dat) >= 5){
      print(dat[1:5, 1:5])
    } else {
      print(head(dat))
    }
  }

  return(dat)
}

linear_kernel <- function(X) {
  # spot x celltype
  std <- scale(X)
  K <- tcrossprod(as.matrix(X)) # K = np.dot(X, X.T) # K <- as.matrix(std) %*% as.matrix(t(std))
  # return(K / max(K))            # return K / K.max()
  return(K)
}

# dist_kernel <- function(X, l) {
#   R2 <- as.matrix(dist(X) ** 2)
#   K <- exp(-R2 / (2 * l ** 2))
#   K
# }

spatial_kernel <- function(normalized_expr, coord, kerneltype = "gaussian", bandwidthtype = "SJ", bandwidth.set.by.user = NULL, sparseKernel = FALSE, sparseKernel_tol = 1e-20, sparseKernel_ncore = 1) {
  # cal spatial Kernel using SpatialPCA r package
  # https://lulushang.org/SpatialPCA_Tutorial/slideseq.html

  # normalized_expr: g x n (we used SCTransform)
  # coord: n x 2 (i.e., x and y)

  # The type of bandwidth to be used in Gaussian kernel,
  #   1. "SJ" for Sheather & Jones (1991) method (usually used in small size datasets),
  #   2. "Silverman" for Silverman's ‘rule of thumb’ method (1986)(usually used in large size datasets).

  # scale expr data to calculate "bandwidth"
  expr <- as.matrix(as.data.frame(t(scale(t(normalized_expr)))))

  if (is.null(bandwidth.set.by.user)) {
    bandwidth <- SpatialPCA::bandwidth_select(expr, method = bandwidthtype)
    cat(paste("## cal bandwidth: ", bandwidth, "\n"), sep = "")
  } else {
    bandwidth <- bandwidth.set.by.user
    cat(paste("## select bandwidth by user: ", bandwidth, "\n"), sep = "")
  }

  # scale coordinate data
  location_normalized <- scale(coord)

  if (sparseKernel == FALSE) {
    cat(paste("## Calculating kernel matrix\n"))
    kernelmat <- SpatialPCA::kernel_build(kerneltype = kerneltype, location = location_normalized, bandwidth = bandwidth)
  } else if (sparseKernel == TRUE) {
    cat(paste("## Calculating sparse kernel matrix\n"))
    kernelmat <- SpatialPCA::kernel_build_sparse(
      kerneltype = kerneltype,
      location = location_normalized, bandwidth = bandwidth,
      tol = sparseKernel_tol, ncores = sparseKernel_ncore
    )
  }
  return(kernelmat)
}
