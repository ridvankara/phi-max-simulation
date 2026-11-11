# PHI-max: Composite P-Hacking Detection Index

This repository accompanies the manuscript:

> Kara, R. (2026). *Design Principles for Composite P-Hacking Detection Indices with Complementary Components: A Monte Carlo Simulation Study.* Submitted to *Journal of Statistical Computation and Simulation*.

## Repository scope

This is a **reproducibility package** for the analyses reported in the manuscript. It contains:

- the labelled simulation grid dataset (n = 44,000 synthetic studies) used in all analyses;
- the R scripts that reproduce the inter-component correlation results (Section 3.4), the broad null size validation (Section 2.9), the benchmark against the CS2B test (Section 3.8), and the AUC / paired-bootstrap / p-curve-ablation analyses (Sections 3.5–3.6, Tables 8 and 12);
- the 11 figures included in the manuscript.

The full simulation pipeline (null-distribution generation, hacking-strategy operationalisation, mid-rank ECDF standardisation, weighted-composite optimisation, and PHI-max evaluation) was implemented in **R 4.5.2** and is described step by step in Sections 2.3–2.10 of the manuscript. The simulation grid `data/sim_grid_results.csv` is the direct output of that pipeline; the downstream analyses below are reproduced from it.

## Repository structure

```
phi-max-simulation/
├── README.md
├── LICENSE                             ← MIT
├── R/
│   ├── correlation_analysis.R          ← Inter-component correlations (Section 3.4; Tables 4–5)
│   ├── size_validation.R               ← Broad null size validation (Section 2.9; Table 10)
│   ├── cs2b_benchmark.R                ← Benchmark vs. CS2B (Section 3.8; Table 11)
│   └── auc_and_ablation.R              ← Pooled/condition-wise AUC, paired bootstrap, p-curve ablation (Sections 3.5–3.6; Tables 8, 12)
├── data/
│   ├── sim_grid_results.csv            ← 44,000 labelled simulated studies
│   └── phi_max_size_table_v2.csv       ← 48-cell size-validation results (Table 10)
├── figures/                            ← 11 final manuscript figures
└── session_info.txt                    ← R version and package versions
```

## Dependencies

Base R plus the **`phack`** package (Kranz, 2024), which wraps the conditional
chi-squared test of Elliott, Kudrin & Wüthrich (2022) used by `cs2b_benchmark.R`:

```r
install.packages("phack",
  repos = c("https://skranz.r-universe.dev", "https://cloud.r-project.org"))
```

## Reproducing the analyses

```r
source("R/correlation_analysis.R")   # Section 3.4
source("R/size_validation.R")        # Section 2.9
source("R/cs2b_benchmark.R")         # Section 3.8 (Table 11)
source("R/auc_and_ablation.R")       # Sections 3.5–3.6, Tables 8 and 12
```

`cs2b_benchmark.R` is self-contained. `auc_and_ablation.R` reads
`data/sim_grid_results.csv`; set the column-name variables at the top of the
script to match the dataset before running.

## Citation

If you use this code or data, please cite the manuscript above.

## License

MIT License — see `LICENSE`.
