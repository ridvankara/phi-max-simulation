# =============================================================================
# PHI-max: AUC comparison, paired-bootstrap test, and p-curve ablation
# Reproduces the pooled AUCs (Sections 3.5-3.6), the paired-bootstrap comparison
# of PHI-max vs. the weighted composite (Section 3.6), the condition-wise AUCs
# (Table 8), and the p-curve component-ablation (Table 12).
#
# Input: data/sim_grid_results.csv  (44,000 labelled studies)
# Set the column names below to match your file, then run top-to-bottom.
# =============================================================================
set.seed(20261109)

## ---- column mapping: EDIT to match sim_grid_results.csv --------------------
COL_COND   <- "condition"   # condition label (clean_null, optstop_high, out_high, ...)
COL_K      <- "k"           # number of p-values per study
COL_HACKED <- "hacked"      # 0/1
COL_PCURVE <- "pcurve"      # RAW p-curve score  (change to your column name)
COL_TIVA   <- "tiva"        # RAW TIVA score
COL_CALIP  <- "caliper"     # RAW Caliper score
## If your file already contains mid-rank-standardised columns, set these and
## skip the standardisation block (USE_MR <- TRUE):
USE_MR <- FALSE
COL_PCURVE_MR <- "pcurve_mr"; COL_TIVA_MR <- "tiva_mr"; COL_CALIP_MR <- "caliper_mr"

d <- read.csv("data/sim_grid_results.csv", stringsAsFactors = FALSE)

## ---- mid-rank ECDF standardisation (from the clean null, per test and k) ---
midrank_std <- function(score, k, calib_by_k) {
  v <- calib_by_k[[as.character(k)]]; n <- length(v)
  (findInterval(score, v, left.open = TRUE) + findInterval(score, v)) / (2 * n)
}
if (USE_MR) {
  d$pcurve_mr <- d[[COL_PCURVE_MR]]; d$tiva_mr <- d[[COL_TIVA_MR]]; d$caliper_mr <- d[[COL_CALIP_MR]]
} else {
  ks <- sort(unique(d[[COL_K]]))
  null_rows <- d[[COL_COND]] == "clean_null"
  for (raw in c(COL_PCURVE, COL_TIVA, COL_CALIP)) {
    calib <- lapply(ks, function(k) sort(d[[raw]][null_rows & d[[COL_K]] == k])); names(calib) <- as.character(ks)
    nm <- c(pcurve = "pcurve_mr", tiva = "tiva_mr", caliper = "caliper_mr")[
            c(COL_PCURVE="pcurve", COL_TIVA="tiva", COL_CALIP="caliper")[raw]]
    d[[nm]] <- mapply(function(s, k) midrank_std(s, k, calib), d[[raw]], d[[COL_K]])
  }
}
d$PHI_max      <- pmax(d$pcurve_mr, d$tiva_mr, d$caliper_mr)
d$PHI_nopcurve <- pmax(d$tiva_mr, d$caliper_mr)           # ablation: drop p-curve

## ---- AUC (Mann-Whitney, tie-corrected) ------------------------------------
auc <- function(score, y) {
  pos <- score[y == 1]; neg <- score[y == 0]
  pos <- pos[is.finite(pos)]; neg <- neg[is.finite(neg)]
  n1 <- as.double(length(pos)); n0 <- as.double(length(neg)); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(c(pos, neg)); (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ---- balanced train/test split + learned composite weights ----------------
hk <- d[d[[COL_HACKED]] == 1, ]; cl <- d[d[[COL_COND]] == "clean_null", ]
set.seed(1); cl <- cl[sample(nrow(cl), nrow(hk)), ]
pool <- rbind(hk, cl); pool$y <- pool[[COL_HACKED]]
tr <- sample(nrow(pool), 0.7 * nrow(pool)); train <- pool[tr, ]; test <- pool[-tr, ]
fit <- glm(y ~ pcurve_mr + tiva_mr + caliper_mr, data = train, family = binomial)
w <- pmax(coef(fit)[c("pcurve_mr", "tiva_mr", "caliper_mr")], 0); w <- w / sum(w)
test$composite <- as.numeric(w[1]*test$pcurve_mr + w[2]*test$tiva_mr + w[3]*test$caliper_mr)

## ---- pooled AUCs + paired-bootstrap delta (Section 3.6) --------------------
A_c <- auc(test$composite, test$y); A_p <- auc(test$PHI_max, test$y)
B <- 2000; delt <- numeric(B); idx <- seq_len(nrow(test))
set.seed(3)
for (b in 1:B) { s <- sample(idx, replace = TRUE); tb <- test[s, ]
  delt[b] <- auc(tb$PHI_max, tb$y) - auc(tb$composite, tb$y) }
ci <- quantile(delt, c(.025, .975)); pval <- mean(delt <= 0)
cat(sprintf("\n[Sec 3.6] Composite AUC=%.3f | PHI-max AUC=%.3f | dAUC=%.3f | 95%% CI [%.3f, %.3f] | one-sided p=%.4f\n",
            A_c, A_p, A_p - A_c, ci[1], ci[2], pval))

## ---- condition-wise AUC (Table 8) -----------------------------------------
neg <- test[test[[COL_COND]] == "clean_null", ]
conds <- setdiff(unique(test[[COL_COND]]), c("clean_null", "clean_mixed"))
tab8 <- do.call(rbind, lapply(conds, function(cc) {
  pos <- test[test[[COL_COND]] == cc, ]
  data.frame(condition = cc,
    composite = round(auc(c(pos$composite, neg$composite), c(rep(1,nrow(pos)), rep(0,nrow(neg)))), 3),
    pcurve    = round(auc(c(pos$pcurve_mr, neg$pcurve_mr), c(rep(1,nrow(pos)), rep(0,nrow(neg)))), 3),
    tiva      = round(auc(c(pos$tiva_mr,   neg$tiva_mr),   c(rep(1,nrow(pos)), rep(0,nrow(neg)))), 3),
    caliper   = round(auc(c(pos$caliper_mr,neg$caliper_mr),c(rep(1,nrow(pos)), rep(0,nrow(neg)))), 3)) }))
cat("\n[Table 8] condition-wise AUC vs clean null:\n"); print(tab8)

## ---- p-curve ablation (Table 12) ------------------------------------------
tab12 <- do.call(rbind, lapply(conds, function(cc) {
  pos <- test[test[[COL_COND]] == cc, ]
  a3 <- auc(c(pos$PHI_max,      neg$PHI_max),      c(rep(1,nrow(pos)), rep(0,nrow(neg))))
  a2 <- auc(c(pos$PHI_nopcurve, neg$PHI_nopcurve), c(rep(1,nrow(pos)), rep(0,nrow(neg))))
  data.frame(condition = cc, PHImax_3comp = round(a3,3), PHImax_no_pcurve = round(a2,3),
             delta_pcurve = round(a3 - a2, 3)) }))
cat("\n[Table 12] p-curve ablation:\n"); print(tab12)

write.csv(tab8,  "data/table8_condition_auc.csv",  row.names = FALSE)
write.csv(tab12, "data/table12_pcurve_ablation.csv", row.names = FALSE)
