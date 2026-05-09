# Mendelian-randomization

This repository contains R code, simulation outputs, figures, and supplementary tables for the simulation study comparing the following Mendelian randomization methods:

* Wald
* 2SPS
* 2SRI
* GMM (corrected)
* IV-MVB (corrected)

The repository accompanies the methodological correction and simulation study described in the submitted Letter to the Editor.

---

# Requirements

R ≥ 4.4

Required R packages:

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "stringr",
  "readr", "gt", "sandwich", "gridExtra",
  "readxl", "purrr", "tibble"
))
```

---

# How to Run

Open R and run:

```r
source("scripts/run_all.R")
```

---

# Repository Structure

```text
.
├── IVMVB_simulation_main.R
├── scripts/
│   └── run_all.R
├── results/
│   ├── figures/
│   └── tables/
```

---

# Simulation Outputs

Simulation figures and supplementary tables used in the manuscript are included in:

```text
results/figures/
results/tables/
```

---

# Reproducibility

* Random seeds are fixed for reproducibility.
* All simulation results can be regenerated using the provided scripts.

---

# Contact

For questions regarding the code or simulations, please contact the authors.
