# Mendelian-randomization

This repository contains R code for the simulation study comparing the following methods:

* Wald
* 2SPS
* 2SRI
* GMM
* IV-MVB

---

## Requirements

R ≥ 4.4

Install required packages:

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "stringr",
  "readr", "gt", "sandwich", "gridExtra",
  "readxl", "purrr", "tibble"
))
```

---

## How to Run

Open R and run:

```r
source("scripts/run_all.R")
```

---

## Notes

* Output (figures and tables) will be saved in:

```
results/figures/
results/tables/
```

---

## Real Data Analysis (Optional)

The CLEAR dataset analysis is included but **disabled by default**.

To enable:

```r
run_clear_analysis <- TRUE
```

Note: The CLEAR dataset is not included in this repository.

---

## Structure

```
.
├── ivmvb_simulation_main.R
├── scripts/
│   └── run_all.R
├── results/
│   ├── figures/
│   └── tables/
```

---

## Reproducibility

* Random seed is fixed for reproducibility.
* All results can be regenerated from the provided scripts.

---

## Contact

For questions, please contact the authors.
