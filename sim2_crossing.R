# sim2_crossing.R
# Study 2: extended simulation for the crossing ratio (CR). See PROTOCOL.md, Study 2.
# Part A: single items (K = 4, 5, 7; normal, logistic, skewed latent; n = 30 to 1,000).
# Part B: summed scales of 1, 3, and 5 items.
# Run by copy and paste into the R console. Output folder is created automatically.
# Required package: ordinal.

library(ordinal)
library(parallel)

out_dir <- "out2"
n_rep   <- 1000
B       <- 999
n_cores <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)
seed0   <- 404
part    <- c("A", "B")         # parts to run
cellsA_run <- NULL             # NULL runs all Part A cells; a range such as 1:1080 runs a chunk
tag     <- ""                  # suffix for chunk output files, e.g. "_1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------- latent distributions (mean 0, SD 1) ----------------
cdf <- list(
  normal   = function(z) pnorm(z),
  logistic = function(z) plogis(z / (sqrt(3) / pi)),
  skewed   = function(z) pchisq(4 + sqrt(8) * z, df = 4))   # standardized chi-square, 4 df

thresholds <- list(
  "4" = list("1" = c(-1, 0, 1),
             "2" = c(-0.5, 0, 0.5),
             "3" = c(-1.6, -1.0, -0.3),
             "4" = c(-1.0, -0.9, 0.8)),
  "5" = list("1" = c(-1.5, -0.5, 0.5, 1.5),
             "2" = c(-1.2, -0.3, 0.3, 1.2),
             "3" = c(-2.0, -1.3, -0.7, -0.2),
             "4" = c(-1.0, -0.2, -0.1, 1.0)),
  "7" = list("1" = c(-1.5, -0.9, -0.3, 0.3, 0.9, 1.5),
             "2" = c(-1.6, -0.9, -0.15, 0.15, 0.9, 1.6),
             "3" = c(-2.2, -1.7, -1.2, -0.8, -0.4, -0.1),
             "4" = c(-1.4, -0.7, -0.05, 0.05, 0.7, 1.4)))

cat_probs <- function(tau, mu, sd, dist) diff(c(0, cdf[[dist]]((tau - mu) / sd), 1))

# ---------------- crossing ratio helpers ----------------
exceed <- function(counts, n) {
  K <- NROW(counts)
  U <- outer(1:(K - 1), 1:K, function(k, j) as.numeric(j > k))
  (U %*% counts) / n
}
cr_cols <- function(cmat, s) {
  pos <- colSums(pmax(s * cmat, 0)); neg <- colSums(pmax(-s * cmat, 0))
  ifelse(neg <= 1e-12, Inf, pos / neg)
}
pop_cr <- function(pA, pB) {
  cd <- exceed(pB, 1) - exceed(pA, 1); D <- sum(cd)
  if (abs(D) < 1e-12) return(c(D = 0, cr = 1))
  c(D = D, cr = cr_cols(matrix(cd), sign(D)))
}
cr_boot <- function(ca, cb) {
  nA <- sum(ca); nB <- sum(cb)
  cd <- exceed(cb, nB) - exceed(ca, nA); D <- sum(cd)
  s <- if (D == 0) 1 else sign(D)
  cr <- if (D == 0) 1 else cr_cols(matrix(cd), s)
  cbm <- exceed(rmultinom(B, nB, cb / nB), nB) - exceed(rmultinom(B, nA, ca / nA), nA)
  crb <- cr_cols(cbm, s)
  c(D = D, cr = cr, lb = unname(quantile(crb, 0.05, type = 1)))
}

welch_counts <- function(ca, cb, x) {
  na <- sum(ca); nb <- sum(cb)
  ma <- sum(x * ca) / na; mb <- sum(x * cb) / nb
  va <- sum(ca * (x - ma)^2) / (na - 1); vb <- sum(cb * (x - mb)^2) / (nb - 1)
  se2 <- va / na + vb / nb
  if (se2 <= 0) return(c(t = NA, p = NA))
  t <- (mb - ma) / sqrt(se2)
  df <- se2^2 / ((va / na)^2 / (na - 1) + (vb / nb)^2 / (nb - 1))
  c(t = t, p = 2 * pt(-abs(t), df))
}

# Wilcoxon-Mann-Whitney from counts (normal approximation with tie correction)
wmw_counts <- function(ca, cb) {
  na <- sum(ca); nb <- sum(cb); N <- na + nb; tot <- ca + cb
  rk <- cumsum(tot) - (tot - 1) / 2                  # mid-ranks per category
  W <- sum(cb * rk) - nb * (nb + 1) / 2
  mu <- na * nb / 2
  v <- na * nb / 12 * ((N + 1) - sum(tot^3 - tot) / (N * (N - 1)))
  if (v <= 0) return(NA)
  z <- (W - mu - 0.5 * sign(W - mu)) / sqrt(v)
  2 * pnorm(-abs(z))
}

# cumulative probit fits on the aggregated table (weights = cell counts)
clm_tests <- function(ca, cb) {
  K <- length(ca); keep <- which(ca + cb > 0); k_obs <- length(keep)
  if (k_obs < 3) return(c(p_clm = NA, p_eq = NA))
  w <- c(ca[keep], cb[keep]); ok <- w > 0
  y <- factor(rep(seq_len(k_obs), 2), ordered = TRUE)[ok]
  X <- cbind("(Intercept)" = 1, g = rep(0:1, each = k_obs))[ok, , drop = FALSE]
  w <- w[ok]
  tot <- (ca + cb)[keep]
  ll0 <- sum(tot * log(tot / sum(tot)))          # null model log-likelihood (closed form)
  f1 <- try(suppressWarnings(clm.fit(y, X, weights = w, link = "probit")), silent = TRUE)
  p_clm <- if (inherits(f1, "try-error")) NA else
    pchisq(max(0, 2 * (f1$logLik - ll0)), 1, lower.tail = FALSE)
  p_eq <- NA
  if (k_obs >= 4 && !inherits(f1, "try-error")) {
    fe <- try(suppressWarnings(clm.fit(y, X, weights = w, link = "probit", threshold = "equidistant")), silent = TRUE)
    if (!inherits(fe, "try-error"))
      p_eq <- pchisq(max(0, 2 * (f1$logLik - fe$logLik)), k_obs - 3, lower.tail = FALSE)
  }
  c(p_clm = p_clm, p_eq = p_eq)
}

# ======================= Part A: single items =======================
designA <- rbind(
  expand.grid(dist = names(cdf), K = c(4, 5, 7), pattern = 1:4, delta = c(0, 0.2, 0.5),
              sdr = c(1, 1.25, 1.5, 2), nA = c(30, 100, 300, 1000), stringsAsFactors = FALSE),
  expand.grid(dist = names(cdf), K = c(4, 5, 7), pattern = 1:4, delta = c(0, 0.2, 0.5),
              sdr = c(1, 1.25, 1.5, 2), nA = 100, stringsAsFactors = FALSE))
designA$nB <- designA$nA
designA$nB[(nrow(designA) - 431):nrow(designA)] <- 300
designA$cell <- seq_len(nrow(designA))

popA <- do.call(rbind, lapply(seq_len(nrow(designA)), function(i) {
  r <- designA[i, ]; tau <- thresholds[[as.character(r$K)]][[as.character(r$pattern)]]
  pA <- cat_probs(tau, 0, 1, r$dist); pB <- cat_probs(tau, r$delta, r$sdr, r$dist)
  pc <- pop_cr(pA, pB)
  data.frame(cell = r$cell, pop_D = pc["D"], pop_cr = pc["cr"],
             latent_disagree = r$delta > 0 & pc["D"] < 0,
             pA = paste(round(pA, 4), collapse = " "), pB = paste(round(pB, 4), collapse = " "))
}))

runA <- function(i) {
  r <- designA[i, ]; K <- r$K
  tau <- thresholds[[as.character(K)]][[as.character(r$pattern)]]
  pA <- cat_probs(tau, 0, 1, r$dist); pB <- cat_probs(tau, r$delta, r$sdr, r$dist)
  set.seed(seed0 * 100000 + r$cell)
  res <- t(replicate(n_rep, {
    ca <- tabulate(sample.int(K, r$nA, TRUE, pA), K)
    cb <- tabulate(sample.int(K, r$nB, TRUE, pB), K)
    c(welch_counts(ca, cb, 1:K), p_w = wmw_counts(ca, cb), clm_tests(ca, cb), cr_boot(ca, cb))
  }))
  data.frame(cell = r$cell, rep = seq_len(n_rep), res)
}

# ======================= Part B: summed scales =======================
lambda <- sqrt(0.6)                      # loading of each item on the common factor
designB <- expand.grid(m = c(1, 3, 5), pattern = c(1, 3), delta = c(0.2, 0.5),
                       sdr = c(1, 1.5, 2), nA = c(100, 300, 1000))
designB$cell <- seq_len(nrow(designB))
tauB <- thresholds[["5"]]

# exact population distribution of the sum score by numerical integration over the factor
sum_probs <- function(m, tau, mu, sd) {
  f <- seq(-8, 8, by = 0.01); wf <- dnorm(f) * 0.01
  out <- 0
  for (j in seq_along(f)) {
    eta <- lambda * (mu + sd * f[j])
    p <- diff(c(0, pnorm((tau - eta) / sqrt(1 - lambda^2)), 1))
    s <- 1; for (k in seq_len(m)) s <- convolve(s, rev(p), type = "open")
    out <- out + wf[j] * s
  }
  out / sum(out)
}
popB <- do.call(rbind, lapply(seq_len(nrow(designB)), function(i) {
  r <- designB[i, ]; tau <- tauB[[as.character(r$pattern)]]
  pA <- sum_probs(r$m, tau, 0, 1); pB <- sum_probs(r$m, tau, r$delta, r$sdr)
  pc <- pop_cr(pA, pB)
  data.frame(cell = r$cell, pop_D = pc["D"], pop_cr = pc["cr"])
}))

runB <- function(i) {
  r <- designB[i, ]; tau <- tauB[[as.character(r$pattern)]]; L <- r$m * 4 + 1
  set.seed(seed0 * 100000 + 50000 + r$cell)
  draw <- function(n, mu, sd) {
    f <- rnorm(n, mu, sd)
    y <- sapply(seq_len(r$m), function(k)
      findInterval(lambda * f + sqrt(1 - lambda^2) * rnorm(n), tau))  # 0..4 per item
    tabulate(rowSums(matrix(y, n)) + 1, L)
  }
  res <- t(replicate(n_rep, {
    ca <- draw(r$nA, 0, 1); cb <- draw(r$nA, r$delta, r$sdr)
    c(welch_counts(ca, cb, 0:(L - 1)), cr_boot(ca, cb))
  }))
  data.frame(cell = r$cell, rep = seq_len(n_rep), res)
}

# ---------------- run ----------------
if ("A" %in% part) {
  write.csv(merge(designA, popA, by = "cell"), file.path(out_dir, "cellsA.csv"), row.names = FALSE)
  t0 <- Sys.time()
  idx <- if (is.null(cellsA_run)) seq_len(nrow(designA)) else cellsA_run
  repsA <- do.call(rbind, mclapply(idx, runA, mc.cores = n_cores))
  write.csv(repsA, gzfile(file.path(out_dir, paste0("repsA", tag, ".csv.gz"))), row.names = FALSE)
  cat("Part A:", nrow(repsA), "samples in", round(difftime(Sys.time(), t0, units = "mins"), 1), "minutes\n")
}
if ("B" %in% part) {
  write.csv(merge(designB, popB, by = "cell"), file.path(out_dir, "cellsB.csv"), row.names = FALSE)
  t0 <- Sys.time()
  repsB <- do.call(rbind, mclapply(seq_len(nrow(designB)), runB, mc.cores = n_cores))
  write.csv(repsB, gzfile(file.path(out_dir, "repsB.csv.gz")), row.names = FALSE)
  cat("Part B:", nrow(repsB), "samples in", round(difftime(Sys.time(), t0, units = "mins"), 1), "minutes\n")
}
