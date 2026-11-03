# analysis2_crossing.R
# Study 2: rates by cell, hypotheses H7 to H13, and summed-scale results. See PROTOCOL.md.
# Run after sim2_crossing.R, by copy and paste into the R console.

out_dir <- "out2"
alpha <- 0.05
cellsA <- read.csv(file.path(out_dir, "cellsA.csv"))
fA <- list.files(out_dir, pattern = "^repsA.*\\.csv\\.gz$", full.names = TRUE)
repsA <- do.call(rbind, lapply(fA, function(f) read.csv(gzfile(f))))
repsA <- merge(repsA, cellsA[, c("cell", "pop_cr")], by = "cell")
cat("Part A samples:", nrow(repsA), "in", length(unique(repsA$cell)), "cells\n")

rate <- function(x) mean(x, na.rm = TRUE)
cnt  <- function(x) sum(x, na.rm = TRUE)
nok  <- function(x) sum(!is.na(x))

cell_rates <- function(d) {
  sig <- !is.na(d$p) & d$p < alpha
  pc <- d$pop_cr[1]; cov_ok <- is.finite(pc)
  data.frame(cell = d$cell[1],
             rej_t = rate(d$p < alpha), rej_w = rate(d$p_w < alpha),
             rej_clm = rate(d$p_clm < alpha), rej_eq = rate(d$p_eq < alpha),
             flag15 = rate(d$lb < 1.5), flag2 = rate(d$lb < 2), flag3 = rate(d$lb < 3),
             cover = if (cov_ok) rate(d$lb <= pc) else NA,
             n_sig = sum(sig), x_sigflag = sum(sig & d$lb < 2),
             sig_flag2 = if (sum(sig) > 0) mean(d$lb[sig] < 2) else NA,
             x_t = cnt(d$p < alpha), n_t = nok(d$p),
             x_w = cnt(d$p_w < alpha), n_w = nok(d$p_w),
             x_eq = cnt(d$p_eq < alpha), n_eq = nok(d$p_eq),
             x_flag2 = cnt(d$lb < 2), x_cover = if (cov_ok) cnt(d$lb <= pc) else NA, n = nrow(d))
}
ratesA <- do.call(rbind, lapply(split(repsA, repsA$cell), cell_rates))
ratesA <- merge(cellsA, ratesA, by = "cell")
ratesA$type <- with(ratesA, ifelse(delta == 0 & sdr == 1, "null",
                           ifelse(is.finite(pop_cr) & pop_cr < 2, "fragile", "safe")))
write.csv(ratesA, file.path(out_dir, "ratesA.csv"), row.names = FALSE)

# ---------------- hypotheses ----------------
p_le <- function(x, n, c) pbinom(x, n, c)
p_ge <- function(x, n, c) pbinom(x - 1, n, c, lower.tail = FALSE)
hyp <- list()
add <- function(id, desc, sub, xcol, ncol, lo = NA, hi = NA) {
  x <- sub[[xcol]]; n <- sub[[ncol]]; keep <- n > 0; x <- x[keep]; n <- n[keep]; cells <- sub$cell[keep]
  if (length(x) == 0) { cat(id, ": no cells available\n"); return(invisible(NULL)) }
  est <- x / n; p <- rep(0, length(x))
  if (!is.na(lo)) p <- pmax(p, p_ge(x, n, lo))
  if (!is.na(hi)) p <- pmax(p, p_le(x, n, hi))
  ok <- (is.na(lo) | est >= lo) & (is.na(hi) | est <= hi)
  w <- if (all(ok)) which.max(p) else which(!ok)[which.max(abs(est[!ok] - ifelse(is.na(lo), hi, lo)))]
  hyp[[length(hyp) + 1]] <<- data.frame(
    hypothesis = id, criterion = desc, cells = length(x), cells_failing = sum(!ok),
    min_rate = round(min(est), 3), max_rate = round(max(est), 3),
    worst_cell = cells[w], worst_rate = round(est[w], 3), worst_n = n[w], p_worst = signif(p[w], 3),
    decision = if (all(ok)) "Supported" else "Not supported")
}
h7 <- subset(ratesA, type == "null" & nA >= 100)
add("H7 (Welch t)", "Type I error 0.025 to 0.075, n >= 100", h7, "x_t", "n_t", lo = 0.025, hi = 0.075)
add("H7 (Wilcoxon)", "Type I error 0.025 to 0.075, n >= 100", h7, "x_w", "n_w", lo = 0.025, hi = 0.075)
add("H8", "CR flag rate <= 0.05 (population CR >= 10, n = 1000)",
    subset(ratesA, type != "null" & (!is.finite(pop_cr) | pop_cr >= 10) & nA == 1000 & nB == 1000), "x_flag2", "n", hi = 0.05)
add("H9", "CR flag rate >= 0.80 (population CR < 1.5, n >= 300)",
    subset(ratesA, type != "null" & pop_cr < 1.5 & nA >= 300 & nA == nB), "x_flag2", "n", lo = 0.80)
add("H10", "Equidistance rejection >= 0.80 (normal, SD ratio 1, delta 0.5, n 1000, patterns 2-4)",
    subset(ratesA, dist == "normal" & sdr == 1 & delta == 0.5 & nA == 1000 & nB == 1000 & pattern %in% 2:4), "x_eq", "n_eq", lo = 0.80)
add("H11", "Equidistance rejection <= 0.50 (fragile, pattern 1, normal, n 1000)",
    subset(ratesA, type == "fragile" & pattern == 1 & dist == "normal" & nA == 1000 & nB == 1000), "x_eq", "n_eq", hi = 0.50)
add("H12", "Coverage >= 0.925 (1 < population CR < Inf, n >= 300)",
    subset(ratesA, is.finite(pop_cr) & pop_cr > 1 + 1e-9 & nA >= 300 & nA == nB), "x_cover", "n", lo = 0.925)
add("H13a", "Significant t results flagged >= 0.80 (population CR < 1.5, n 1000)",
    subset(ratesA, type != "null" & pop_cr < 1.5 & nA == 1000 & nB == 1000 & n_sig > 0), "x_sigflag", "n_sig", lo = 0.80)
add("H13b", "Significant t results flagged <= 0.20 (population CR >= 10, n 1000)",
    subset(ratesA, type != "null" & (!is.finite(pop_cr) | pop_cr >= 10) & nA == 1000 & nB == 1000 & n_sig > 0), "x_sigflag", "n_sig", hi = 0.20)
hyp <- do.call(rbind, hyp)
write.csv(hyp, file.path(out_dir, "hypotheses2.csv"), row.names = FALSE)
print(hyp[, -2], row.names = FALSE)

# ---------------- descriptive summaries ----------------
ratesA$crband <- cut(ifelse(is.finite(ratesA$pop_cr), ratesA$pop_cr, 1e6),
                     c(0, 1.5, 2, 3, 5, 10, Inf), right = FALSE,
                     labels = c("<1.5", "1.5-2", "2-3", "3-5", "5-10", ">=10"))
eqn <- subset(ratesA, type != "null" & nA == nB)
tabF <- aggregate(cbind(flag15, flag2, flag3, rej_eq) ~ crband + nA, data = eqn, FUN = mean)
write.csv(tabF, file.path(out_dir, "flag_by_band.csv"), row.names = FALSE)
cat("\nFlag and equidistance rates by population CR band and n:\n"); print(tabF, digits = 2)
tabS <- aggregate(sig_flag2 ~ crband + nA, data = subset(eqn, n_sig >= 20), FUN = mean)
write.csv(tabS, file.path(out_dir, "sigflag_by_band.csv"), row.names = FALSE)
cat("\nLatent mean and category result disagree (cells, by distribution):\n")
print(table(cellsA$dist[cellsA$latent_disagree]))
dis <- subset(ratesA, latent_disagree & nA == 1000 & nB == 1000)
cat("At n = 1000, Welch t rejection in those cells: median", round(median(dis$rej_t), 3),
    "max", round(max(dis$rej_t), 3), "\n")

# ---------------- Part B ----------------
if (file.exists(file.path(out_dir, "repsB.csv.gz"))) {
  cellsB <- read.csv(file.path(out_dir, "cellsB.csv"))
  repsB <- read.csv(gzfile(file.path(out_dir, "repsB.csv.gz")))
  rB <- do.call(rbind, lapply(split(repsB, repsB$cell), function(d) data.frame(
    cell = d$cell[1], rej_t = rate(d$p < alpha), flag2 = rate(d$lb < 2))))
  rB <- merge(cellsB, rB, by = "cell")
  write.csv(rB, file.path(out_dir, "ratesB.csv"), row.names = FALSE)
  cat("\nPart B population CR by number of items:\n")
  print(reshape(unique(rB[, c("m", "pattern", "delta", "sdr", "pop_cr")]),
                idvar = c("pattern", "delta", "sdr"), timevar = "m", direction = "wide"), digits = 3)
  cat("\nPart B flag rate at n = 1000 by m:\n")
  print(aggregate(flag2 ~ m + pattern + sdr, data = subset(rB, nA == 1000), FUN = mean), digits = 2)
}

# ---------------- reading rule outcomes and threshold sensitivity (paper Table 4 and Section 4.2) ----------------
rr <- merge(repsA[, c("cell", "cr", "lb", "pop_cr")], cellsA[, c("cell", "nA", "nB", "delta", "sdr")], by = "cell")
rr <- subset(rr, nA == nB & !(delta == 0 & sdr == 1))
rr$band <- cut(ifelse(is.finite(rr$pop_cr), rr$pop_cr, 1e6), c(0, 1.5, 2, 3, 5, 10, Inf), right = FALSE,
               labels = c("<1.5", "1.5-2", "2-3", "3-5", "5-10", ">=10"))
rr$robust  <- rr$lb >= 2
rr$depends <- !rr$robust & rr$cr < 2 & rr$nA >= 300
rule <- aggregate(cbind(robust, depends) ~ band + nA, data = rr, FUN = mean)
rule$undecided <- 1 - rule$robust - rule$depends
write.csv(rule, file.path(out_dir, "rule_outcomes.csv"), row.names = FALSE)
sens <- do.call(rbind, lapply(c(1.5, 2, 3, 5), function(r0) do.call(rbind, lapply(c(30, 100, 300, 1000), function(n) {
  d <- rr[rr$nA == n, ]
  data.frame(r0 = r0, n = n,
             detect = mean(d$lb[is.finite(d$pop_cr) & d$pop_cr < r0] < r0),
             false_flag_no_crossing = mean(d$lb[!is.finite(d$pop_cr)] < r0))
}))))
write.csv(sens, file.path(out_dir, "threshold_sensitivity.csv"), row.names = FALSE)
print(rule, digits = 2); print(sens, digits = 2)
