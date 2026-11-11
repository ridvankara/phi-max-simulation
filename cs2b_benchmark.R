# =============================================================================
# PHI-max: Benchmark against the CS2B test (Elliott, Kudrin & Wüthrich, 2022)
# Reproduces Table 11 and the small-corpus size-breakdown reported in Section 3.8.
#
# CS2B = conditional chi-squared test of 2-monotonicity and bounds, computed with
#        the authors' own implementation via the `phack` package (Kranz, 2024):
#        phack_test_cox_shi(p, K = 2, use_bounds = TRUE, J = 15, p_min = 0, p_max = 0.15)
#
# Requires: install.packages("phack",
#             repos = c("https://skranz.r-universe.dev", "https://cloud.r-project.org"))
# =============================================================================
library(phack)
set.seed(20261109)

N_GROUP <- 30           # per-group sample size of the base two-sample z-test
NCP     <- qnorm(0.975) # noncentrality for ~50% power (matches main pipeline)
F_TRUE  <- 0.40         # proportion of true effects in the heterogeneous condition
MB      <- 2L           # Caliper minimum band count (matches main pipeline)

# ---- CS2B decision (returns 1 = reject "no p-hacking" at 0.05, else 0) -------
cs2b_rej <- function(p) {
  v <- tryCatch(phack_test_cox_shi(p, K = 2, use_bounds = TRUE,
                                   J = 15, p_min = 0, p_max = 0.15),
                error = function(e) NA)
  if (is.null(v) || is.na(v) || v > 1) return(0L)   # 999/888 = non-convergence -> non-rejection
  as.integer(v < 0.05)
}

# ---- Caliper decision (corpus-level comparator) -----------------------------
caliper_rej <- function(p) {
  lo <- sum(p >= 0.040 & p < 0.045); up <- sum(p >= 0.045 & p < 0.050); n <- lo + up
  if (n < MB) return(0L)
  as.integer(binom.test(up, n, 0.5, alternative = "greater")$p.value < 0.05)
}

# ---- generators: one reported p-value per researcher ------------------------
opt_stop_p <- function(max_extra, delta = 0, n_start = 20, n_step = 5) {
  nmax <- n_start + max_extra
  x1 <- rnorm(nmax); x2 <- rnorm(nmax) + delta
  cs1 <- cumsum(x1); cs2 <- cumsum(x2)
  looks <- seq(n_start, nmax, by = n_step)
  z <- (cs2[looks]/looks - cs1[looks]/looks) / sqrt(2/looks)
  pp <- 2 * (1 - pnorm(abs(z))); hit <- which(pp < 0.05)
  if (length(hit)) pp[hit[1]] else pp[length(pp)]
}
gen <- list(
  clean_null  = function(Np) runif(Np),                                              # homogeneous global null (delta = 0)
  clean_mixed = function(Np) { z <- rnorm(Np) + (runif(Np) < F_TRUE) * NCP; 2*(1-pnorm(abs(z))) }, # heterogeneous non-hacking (delta > 0)
  optstop_high = function(Np) vapply(seq_len(Np), function(i) opt_stop_p(200), 0),   # optional stopping
  minp_high    = function(Np) apply(matrix(runif(Np*10), Np, 10), 1, min)            # minimum approach (specification search)
)

# ---- benchmark across conditions and corpus sizes ---------------------------
R_LIT <- 300
NP    <- c(250, 1000, 3000, 6000)   # corpus sizes; Table 11 reports Np = 3000, 6000
res <- data.frame()
t0 <- Sys.time()
for (Np in NP) for (cc in names(gen)) {
  cs <- ca <- numeric(R_LIT)
  for (r in seq_len(R_LIT)) { p <- gen[[cc]](Np); cs[r] <- cs2b_rej(p); ca[r] <- caliper_rej(p) }
  res <- rbind(res, data.frame(Np = Np, condition = cc,
                               CS2B = round(mean(cs), 3), Caliper = round(mean(ca), 3)))
  cat(sprintf("Np=%5d  %-12s  CS2B=%.3f  Caliper=%.3f  | %.1f min\n",
              Np, cc, mean(cs), mean(ca), as.numeric(Sys.time()-t0, units = "mins")))
}
write.csv(res, "data/cs2b_benchmark_results.csv", row.names = FALSE)
cat("\n=== CS2B benchmark (Table 11 = Np >= 3000 rows) ===\n"); print(res)

# Notes:
#  - At Np >= 3000 CS2B is approximately size-controlling (null ~0.08-0.09; delta>0 ~0.02-0.04)
#    and attains full power against both hacking strategies.
#  - At Np = 250-1000 CS2B over-rejects under the null (too few p-values in [0,0.15]),
#    documenting its large-corpus requirement (Section 3.8).
