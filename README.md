# Rabt — Aggregate Boosting Trees for Ecology & Microbiome

R package for Adaptive Boosting Tree (ABT) modeling based on Gradient Boosting Machine (GBM), with cross-validation for optimal tree selection. Designed for ecology, microbiome, and environmental-community relationship analysis.

---

## Installation

```r
remotes::install_github("Wang-lang/Rabt")
```

## Quick Start

```r
library(Rabt)

data(softcorals)
fit <- ABT(
  formula  = Richness ~ Across + Along + Visibility + Slope + Flow + Wave + Sediment,
  data     = softcorals,
  n.trees  = 100,
  cv.folds = 5
)
abt_plot(fit)
```

---

## Core Functions

### `ABT()` — Model Fitting

Adaptive boosting tree model with k-fold cross-validation to prevent overfitting.

```r
ABT(
  formula, data,
  distribution      = "gaussian",
  n.trees           = 200,
  interaction.depth = 3,
  n.minobsinnode    = NULL,
  shrinkage         = 0.05,
  bag.fraction      = 0.5,
  stratify          = TRUE,
  monitor           = TRUE,
  use               = "best",
  seed              = 0,
  cv.folds          = 5,
  ...
)
```

| Parameter | Description |
|-----------|-------------|
| `distribution` | `"gaussian"`, `"bernoulli"`, `"poisson"`, `"adaboost"`, `"laplace"`, `"coxph"` |
| `n.trees` | Total boosting iterations |
| `interaction.depth` | Tree depth; 1 = no interactions |
| `shrinkage` | Learning rate |
| `bag.fraction` | Subsampling fraction per tree |
| `cv.folds` | Cross-validation folds |
| `stratify` | Stratified CV sampling |

**Returns** an `ABT` object with:
- `gbm.model` — final GBM model
- `best.n.trees` — optimal tree count
- `best.cv.error` — CV error at optimum
- `cv.error` / `cv.error.matrix` — CV error trajectories
- `var.importance` — variable importance table

---

### `micro_data()` — Microbiome Data Preparation

Prepares pairwise difference data from environmental and OTU/ASV tables for ABT analysis.

```r
micro_data(
  env_data,
  otu_data,
  datatype      = c("otu/asv", "matrix"),
  method        = NULL,
  include_pairs = TRUE
)
```

| Parameter | Description |
|-----------|-------------|
| `env_data` | Environmental variables (rows = samples) |
| `otu_data` | OTU/ASV table (rows = taxa, cols = samples) or distance matrix |
| `datatype` | `"otu/asv"` (compute distance) or `"matrix"` (pre-computed) |
| `method` | Distance method for `vegdist`: `"bray"`, `"jaccard"`, `"euclidean"`, `"manhattan"`, `"canberra"`, `"horn"`, `"gower"` |

**Returns** a data frame with `N²` rows (N = samples):
- `Sample_i`, `Sample_j` — pair identifiers
- `env_diff` — environmental differences
- `otu_value` — dissimilarity value

---

### `abt_plot()` — Variable Importance Plot

Publication-ready horizontal bar plot of variable relative importance.

```r
abt_plot(
  fit,
  dict_label = NULL,
  out_plot   = ".",
  plot_title = "Relative Influence of Variables",
  x_label    = "Relative Influence",
  y_label    = "Variable",
  width      = 10, height = 8, dpi = 300
)
```

**Output**: prints to device, saves `ABT.png` and `ABT.pdf`.

---

## Datasets

| Dataset | Description |
|---------|-------------|
| `softcorals` | Soft coral richness and environmental variables |
| `env` | Environmental variables for microbiome analysis |
| `otu` | OTU/ASV count table (taxa × samples) |

---

## Dependencies

`gbm`, `ggplot2`

---

Citation
If you use this package in your research, please cite:
De'ath G. Boosted trees for ecological modeling and prediction. Ecology. 2007 Jan;88(1):243-51.

---

## License

MIT
