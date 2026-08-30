# PHI-max: Composite P-Hacking Detection Index

This repository accompanies the manuscript:

> Kara, R. (2026). *Design Principles for Composite P-Hacking Detection Indices with Complementary Components: A Monte Carlo Simulation Study.* Submitted to *Journal of Statistical Computation and Simulation*.

## Repository scope

This is a **minimal reproducibility package**. It contains:
- The labelled simulation grid dataset (n = 44,000 synthetic studies) used in all analyses reported in the manuscript;
- The correlation analysis script that produces the inter-component correlation results (Section 3.4) — the principal post-review revision of the manuscript;
- The 11 figures included in the manuscript.

The full simulation pipeline (null-distribution generation, hacking-strategy operationalisation, mid-rank ECDF standardisation, weighted-composite optimisation, and PHI-max evaluation) was implemented in R 4.5.2 and is described step by step in Sections 2.3–2.8 of the manuscript. The simulation grid `sim_grid_results.csv` provided here is the direct output of that pipeline; downstream analyses can be reproduced from this file. The full simulation code is available from the corresponding author upon reasonable request.

## Repository structure

```
phi-max-simulation/
├── README.md                         ← This file
├── LICENSE                           ← MIT
├── R/
│   └── correlation_analysis.R        ← Reproduces Section 3.4 correlation tables
├── data/
│   └── sim_grid_results.csv          ← 44,000 simulated studies (labelled)
├── figures/                          ← 11 final figures from the manuscript
└── session_info.txt                  ← R version and package versions used
```

## Dataset description (`data/sim_grid_results.csv`)

44,000 rows × 15 columns. One row per simulated study.

| Column           | Description                                                                |
|------------------|----------------------------------------------------------------------------|
| `condition`      | One of 11 labels (`clean_null`, `clean_mixed`, `optstop_low/med/high`, `sel_low/med/high`, `out_low/med/high`) |
| `label`          | Human-readable label                                                       |
| `strategy`       | One of: `clean`, `optstop`, `selective`, `outcome`                         |
| `hacked`         | Binary indicator (0 = clean control, 1 = hacked)                           |
| `k`              | Number of p-values per study (10, 20, 50, 100)                             |
| `rep`            | Replicate index (1–1000) within condition × k                              |
| `pcurve`         | Raw p-curve score (Stouffer left-skew z)                                   |
| `tiva`           | Raw TIVA score (−log₁₀ p)                                                  |
| `caliper`        | Raw Caliper score (−log₁₀ p)                                               |
| `pcurve_std`     | Standard ECDF-standardised p-curve score (∈ [0,1])                         |
| `tiva_std`       | Standard ECDF-standardised TIVA score                                      |
| `caliper_std`    | Standard ECDF-standardised Caliper score                                   |
| `pcurve_mr`      | Mid-rank ECDF-standardised p-curve score (recommended; see Section 3.3)    |
| `tiva_mr`        | Mid-rank ECDF-standardised TIVA score                                      |
| `caliper_mr`     | Mid-rank ECDF-standardised Caliper score                                   |

The mid-rank ECDF standardisation (`*_mr` columns) is the methodological correction documented in Sections 2.6 and 3.3 of the manuscript and is the recommended input for downstream composite analyses, including PHI-max.

## How to reproduce Section 3.4 (correlation analysis)

1. Clone or download this repository.
2. Open `R/correlation_analysis.R` in R or RStudio.
3. Set the working directory to the repository root, or adjust the data path at the top of the script.
4. Run the script. It produces, in order:
   - **Table 4 Panel A** — pooled inter-component Pearson and Spearman correlation matrices (n = 44,000);
   - **Table 4 Panel B** — within-strategy-family correlation matrices (Clean, Opt.Stop, Min-p Sel+Out);
   - **Table 5** — within-condition correlation table (11 conditions × 3 component pairs);
   - **Section 3.4.3 verification** — selective–outcome operational equivalence check (|Δr| ≤ 0.030);
   - Output CSV files in the working directory.

Run time: < 5 seconds on a standard laptop.

## How to reproduce other manuscript results

The manuscript reports four sets of analyses derived from `sim_grid_results.csv`:

- **Section 3.1 — Null behaviour of component tests** (Figure 1, Table 2): histograms and medians by `k` for `pcurve`, `tiva`, `caliper` filtered to `strategy == "clean"`.
- **Section 3.2 — Component signatures** (Figure 4, Table 3): means of raw `pcurve`, `tiva`, `caliper` grouped by `condition`.
- **Sections 3.3 and 3.4 — Mid-rank correction and complementarity** (Figures 5–6, Table 4–6): comparison of `*_std` vs `*_mr` columns; correlations on `*_mr` columns (this is what `correlation_analysis.R` covers).
- **Sections 3.5–3.7 — Composite optimisation and PHI-max** (Figures 7–11, Tables 7–9): logistic regression and AUC-maximisation on `*_mr` columns; `pmax(pcurve_mr, tiva_mr, caliper_mr)` for PHI-max.

Independent re-implementation from the dataset is straightforward; the manuscript provides full mathematical specifications for each step.

## Software environment

R 4.5.2 (2025-10-31 ucrt) on Windows 11 (64-bit). Principal packages:

- `pROC` (1.19.0.1) — ROC curve construction and AUC computation
- `ggplot2` (4.0.2) — figures
- `RColorBrewer` (1.1-3) — figure palettes
- `dplyr` (1.2.0) — data manipulation

Full session details are in `session_info.txt`.

The correlation analysis (`R/correlation_analysis.R`) requires only base R and runs on R ≥ 4.0.
## Broad Null Size Validation (Revised Manuscript, Section 2.9)
`R/size_validation.R` reproduces the 48-condition broad null validation
reported in Table 10 of the revised manuscript.
Runtime: ~30–60 min (single core). R >= 4.0, base packages only.
Output file: `data/phi_max_size_table_v2.csv`
## Citation

If you use this dataset or analysis in your work, please cite the manuscript above.

## Contact

Rıdvan Kara — Hakkari University, Yüksekova Vocational School, Department of Plant and Animal Production, Hakkari, Türkiye.
ORCID: 0000-0001-6977-4766
Email: ridvankara@hakkari.edu.tr

## Licence

MIT — see `LICENSE`.
