x <- data.frame(matrix(runif(100 * 1e4), ncol = 100))
y <- as.list(x)
medians <- vapply(x, median, numeric(1))

for(i in seq_along(medians)) {
  x[, i] <- x[, i] - medians[i]
}

for(i in 1:5) {
  y[[i]] <- y[[i]] - medians[i]
}