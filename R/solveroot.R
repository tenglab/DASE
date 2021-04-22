solveroot <- function (x, y, y0 = 0) {
  if (is.unsorted(x)) {
    ind <- order(x)
    x <- x[ind]; y <- y[ind]
  }
  z <- y - y0
  # find which part cross zero
  k <- which(z[-1] * z[-length(z)] < 0)
  # find root
  xk <- x[k] - z[k] * (x[k + 1] - x[k]) / (z[k + 1] - z[k])
  xk
}
