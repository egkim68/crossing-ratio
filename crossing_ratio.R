# crossing_ratio.R
# Crossing ratio (CR) for comparing two groups on one Likert item.
# Copy and paste into the R console, then call crossing_ratio(y, group).
#
# y      : answers coded as ordered numbers (for example 1 to 5); NA allowed
# group  : two-level grouping variable (the second level is group B)
# B      : bootstrap resamples for the lower bound
# Lower bound: nonparametric bootstrap within groups (respondents resampled with
#   replacement, drawn as multinomial category counts); approximate one-sided 95% bound.
# r0     : threshold for the flag (default 2)
#
# Reading (from the accompanying paper):
#   lower bound >= r0                        -> direction holds for interval ratios below r0
#   point estimate < r0 and both groups >= 300 -> potentially sensitive to unequal intervals
#   otherwise                                -> undecided

crossing_ratio <- function(y, group, B = 999, r0 = 2, seed = NULL) {
  ok <- !is.na(y) & !is.na(group)
  y <- y[ok]; group <- factor(group[ok])
  if (nlevels(group) != 2) stop("group needs exactly two levels")
  lev <- sort(unique(y)); K <- length(lev)
  if (K < 3) stop("y needs at least three observed categories")
  yi <- match(y, lev)
  ca <- tabulate(yi[group == levels(group)[1]], K)
  cb <- tabulate(yi[group == levels(group)[2]], K)
  nA <- sum(ca); nB <- sum(cb)
  exceed <- function(counts, n) {
    U <- outer(1:(K - 1), 1:K, function(k, j) as.numeric(j > k))
    (U %*% counts) / n
  }
  ratio <- function(cmat, s) {
    pos <- colSums(pmax(s * cmat, 0)); neg <- colSums(pmax(-s * cmat, 0))
    ifelse(neg <= 1e-12, Inf, pos / neg)
  }
  cdiff <- exceed(cb, nB) - exceed(ca, nA)
  D <- sum(cdiff); s <- if (D == 0) 1 else sign(D)
  cr <- if (D == 0) 1 else ratio(matrix(cdiff), s)
  if (!is.null(seed)) set.seed(seed)
  boot <- ratio(exceed(rmultinom(B, nB, cb / nB), nB) - exceed(rmultinom(B, nA, ca / nA), nA), s)
  lb <- unname(quantile(boot, 0.05, type = 1))
  reading <- if (lb >= r0) paste("direction holds for interval ratios below", r0) else
    if (cr < r0 && min(nA, nB) >= 300) "potentially sensitive to unequal intervals" else "undecided"
  out <- list(groups = levels(group), n = c(A = nA, B = nB),
              counts = rbind(A = ca, B = cb), categories = lev,
              cut_point_differences = setNames(round(c(cdiff), 4), paste0("above_", lev[-K])),
              mean_difference = D, crossing_ratio = cr, lower_bound = lb, reading = reading)
  class(out) <- "crossing_ratio"
  out
}

print.crossing_ratio <- function(x, ...) {
  cat("Groups:", x$groups[1], "(A) vs", x$groups[2], "(B); n =", x$n, "\n")
  cat("Cut-point differences (B minus A):", x$cut_point_differences, "\n")
  cat("Mean difference (B minus A):", round(x$mean_difference, 3), "\n")
  cat("Crossing ratio:", round(x$crossing_ratio, 2), " lower 95% bound:", round(x$lower_bound, 2), "\n")
  cat("Reading:", x$reading, "\n")
  invisible(x)
}

# Example (the constructed example in the paper):
# y <- c(rep(3, 2), rep(4, 8), 1, rep(2, 5), 3, rep(5, 13))
# g <- c(rep("men", 10), rep("women", 20))
# crossing_ratio(y, g, seed = 1)
