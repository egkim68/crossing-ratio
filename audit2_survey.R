# audit2_survey.R
# Extended real-data audit of the crossing ratio (Study 2): GSS 2022 and 2024, bfi, and spi.
# Run by copy and paste into the R console. Output folder is created automatically.
# Packages: ordinal, psychTools, and gssr for the GSS. Install gssr with
#   install.packages("gssr", repos = c("https://kjhealy.r-universe.dev", "https://cloud.r-project.org"))

library(ordinal)
library(psychTools)

out_dir <- "out2"
B <- 999
set.seed(505)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gss_rda <- "gss_all.rda"   # used only if gssr is not installed
if (requireNamespace("gssr", quietly = TRUE)) {
  data("gss_all", package = "gssr"); gss <- gss_all
} else {
  e <- new.env(); load(gss_rda, envir = e); gss <- get(ls(e)[1], e)
}

exceed <- function(counts, n) {
  K <- NROW(counts)
  U <- outer(1:(K - 1), 1:K, function(k, j) as.numeric(j > k))
  (U %*% counts) / n
}
cr_cols <- function(cmat, s) {
  pos <- colSums(pmax(s * cmat, 0)); neg <- colSums(pmax(-s * cmat, 0))
  ifelse(neg <= 1e-12, Inf, pos / neg)
}

audit_one <- function(y, g, w = NULL) {   # y: codes; g: 0 = group A, 1 = group B; w: optional weights
  ok <- !is.na(y) & !is.na(g); if (!is.null(w)) ok <- ok & !is.na(w)
  y <- y[ok]; g <- g[ok]; if (!is.null(w)) w <- w[ok]
  lev <- sort(unique(y)); yi <- match(y, lev); K <- length(lev)
  ca <- tabulate(yi[g == 0], K); cb <- tabulate(yi[g == 1], K); nA <- sum(ca); nB <- sum(cb)
  tt <- t.test(y[g == 1], y[g == 0])
  cd <- exceed(cb, nB) - exceed(ca, nA); D <- sum(cd); s <- if (D == 0) 1 else sign(D)
  cr <- if (D == 0) 1 else cr_cols(matrix(cd), s)
  cbm <- exceed(rmultinom(B, nB, cb / nB), nB) - exceed(rmultinom(B, nA, ca / nA), nA)
  lb <- unname(quantile(cr_cols(cbm, s), 0.05, type = 1))
  # equidistance test on the aggregated table
  yy <- factor(rep(1:K, 2), ordered = TRUE); X <- cbind("(Intercept)" = 1, g = rep(0:1, each = K))
  wt <- c(ca, cb); k <- wt > 0
  f1 <- suppressWarnings(clm.fit(yy[k], X[k, ], weights = wt[k], link = "probit"))
  fe <- suppressWarnings(clm.fit(yy[k], X[k, ], weights = wt[k], link = "probit", threshold = "equidistant"))
  p_eq <- pchisq(max(0, 2 * (f1$logLik - fe$logLik)), K - 3, lower.tail = FALSE)
  # weighted CR (point estimate), if weights are given
  cr_w <- NA
  if (!is.null(w)) {
    wa <- sapply(1:K, function(j) sum(w[g == 0 & yi == j])); wb <- sapply(1:K, function(j) sum(w[g == 1 & yi == j]))
    cdw <- exceed(wb, sum(wb)) - exceed(wa, sum(wa)); Dw <- sum(cdw)
    cr_w <- if (Dw == 0) 1 else cr_cols(matrix(cdw), sign(Dw))
  }
  data.frame(K = K, nA = nA, nB = nB, meanA = mean(y[g == 0]), meanB = mean(y[g == 1]),
             t = unname(tt$statistic), p_t = tt$p.value, D = D, cr = cr, lb = lb, p_eq = p_eq,
             cr_weighted = cr_w, exceed_diff = paste(round(cd, 4), collapse = " "))
}

run_source <- function(src, dat, items, groups, w = NULL) {
  do.call(rbind, lapply(items, function(v) {
    if (all(is.na(dat[[v]]))) return(NULL)
    do.call(rbind, lapply(names(groups), function(gn)
      cbind(source = src, item = v, grouping = gn,
            audit_one(as.numeric(dat[[v]]), groups[[gn]], w))))
  }))
}

# ---------------- GSS 2022 and 2024 ----------------
gss_items <- c("satjob", "health", "fefam", "fepresch", "fechld", "fehire", "relpersn",
               "sprtprsn", "finrela", "eqwlth", "helppoor", "helpnot", "helpsick",
               "advfront", "nextgen", "toofast", "meovrwrk", "pillok", "spanking",
               "xmarsex", "homosex")
gss_run <- function(yr) {
  d <- gss[gss$year == yr, ]
  groups <- list(gender = ifelse(as.numeric(d$sex) == 2, 1, 0),          # B = female
                 age    = ifelse(as.numeric(d$age) >= 50, 1, 0),         # B = 50 and over
                 degree = ifelse(as.numeric(d$degree) >= 3, 1, 0))       # B = bachelor or higher
  run_source(paste("GSS", yr), d, gss_items, groups, w = as.numeric(d$wtssps))
}

# ---------------- bfi and spi (psychTools), adults only ----------------
data(bfi, package = "psychTools"); data(spi, package = "psychTools")
bfi <- bfi[!is.na(bfi$age) & bfi$age >= 18, ]
spi <- spi[!is.na(spi$age) & spi$age >= 18, ]
bfi_groups <- list(gender = ifelse(bfi$gender == 2, 1, 0),
                   age    = ifelse(bfi$age >= 30, 1, 0),
                   degree = ifelse(bfi$education >= 4, 1, 0))
spi_groups <- list(gender = ifelse(spi$sex == 2, 1, 0),
                   age    = ifelse(spi$age >= 30, 1, 0),
                   degree = ifelse(spi$education >= 6, 1, 0))
spi_items <- grep("^q_", names(spi), value = TRUE)

audit <- rbind(gss_run(2022), gss_run(2024),
               run_source("bfi", bfi, names(bfi)[1:25], bfi_groups),
               run_source("spi", spi, spi_items, spi_groups))
audit$sig_t <- audit$p_t < 0.05
audit$flag2 <- audit$lb < 2
audit$crossing <- is.finite(audit$cr)
write.csv(audit, file.path(out_dir, "audit2.csv"), row.names = FALSE)

summ <- do.call(rbind, lapply(split(audit, audit$source), function(d) data.frame(
  source = d$source[1], comparisons = nrow(d), significant_t = sum(d$sig_t),
  eq_reject = sum(d$p_eq < 0.05), patterns_cross = sum(d$crossing),
  sig_cross = sum(d$sig_t & d$crossing), sig_flagged = sum(d$sig_t & d$flag2),
  share_sig_flagged = round(mean(d$flag2[d$sig_t]), 3))))
write.csv(summ, file.path(out_dir, "audit2_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)

# replication: GSS 2022 versus 2024, same item and grouping
g22 <- subset(audit, source == "GSS 2022"); g24 <- subset(audit, source == "GSS 2024")
rep_tab <- merge(g22, g24, by = c("item", "grouping"), suffixes = c("_2022", "_2024"))
rep_tab <- rep_tab[, c("item", "grouping", "p_t_2022", "p_t_2024", "cr_2022", "cr_2024",
                       "lb_2022", "lb_2024", "flag2_2022", "flag2_2024")]
write.csv(rep_tab, file.path(out_dir, "audit2_gss_replication.csv"), row.names = FALSE)
cat("\nGSS 2022 vs 2024 comparisons:", nrow(rep_tab),
    " flag agreement:", round(mean(rep_tab$flag2_2022 == rep_tab$flag2_2024), 3),
    " both significant and both unflagged:",
    sum(rep_tab$p_t_2022 < 0.05 & rep_tab$p_t_2024 < 0.05 & !rep_tab$flag2_2022 & !rep_tab$flag2_2024), "\n")

# significant comparisons whose answer patterns cross
sc <- audit[audit$sig_t & audit$crossing, c("source", "item", "grouping", "K", "nA", "nB",
                                            "meanA", "meanB", "p_t", "cr", "lb", "cr_weighted")]
write.csv(sc[order(sc$cr), ], file.path(out_dir, "audit2_sig_crossing.csv"), row.names = FALSE)
cat("Significant comparisons with crossing patterns:", nrow(sc), "\n")
