# Tetraploid Potato GWAS — Mapping Complex Trait Loci

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![GWASpoly](https://img.shields.io/badge/GWASpoly-Polyploid%20GWAS-green.svg)](https://github.com/jendelman/GWASpoly)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

A comprehensive GWAS (Genome-Wide Association Study) pipeline for identifying genetic loci associated with **yield** and **canopy architecture** in a panel of ~282 tetraploid potato (*Solanum tuberosum*) varieties. This project integrates SNP genotype data with multi-year field trial phenotypes and UAV-derived canopy measurements, applying three independent population structure correction approaches to ensure robust marker-trait associations.
Click here for live demo
https://019ce3c9-6a26-c289-c7c2-cb47bf2977f3.share.connect.posit.cloud/
### Highlights

- Tetraploid-aware association mapping with **GWASpoly** (6 genetic models)
- Population structure assessed via **PCA**, **DAPC**, and **STRUCTURE**
- Phenotype estimation using **BLUEs** from mixed-effects models
- Analysis of 5 agronomic/architectural traits across 3 field seasons
- **Interactive R Shiny dashboard** for exploring all results

---

## Repository Layout

```
├── R/
│   ├── 01_phenotype_processing.Rmd   # BLUEs calculation & trait exploration
│   ├── 02_population_structure.Rmd   # PCA, DAPC, STRUCTURE analysis
│   ├── 03_association_mapping.Rmd    # GWASpoly GWAS for all trait × method combinations
│   └── app.R                         # Interactive Shiny dashboard
│
├── data/                             # Input datasets (see data/README.md)
├── docs/                             # Supplementary documentation
├── figures/                          # Generated plots
├── output/                           # Analysis output CSVs
│
├── README.md
├── .gitignore
└── LICENSE
```

---

## Analytical Workflow

### Phase 1 — Phenotype Processing (`01_phenotype_processing.Rmd`)
- Phenotype data for 282 varieties across 2016-2018, two farming systems (Organic & Conventional), 2 replicates
- Missing value imputation using within-variety replicate means
- BLUEs computed via: `Trait ~ Genotype + (1|System) + (1|Year)`
- Model validation through likelihood ratio tests, AIC, and BIC comparisons
- Traits studied: Yield (kg), Max Height (m), Avg Height (m), Canopy Area (m²), Canopy Volume (m³)

### Phase 2 — Population Structure (`02_population_structure.Rmd`)
- SNP filtering: remove markers with ≥20% missing data; PIC > 0.4 for STRUCTURE compatibility
- **PCA**: glPca from adegenet; scree plot elbow at 4 PCs, 5 retained
- **DAPC**: K-means (K=4 via BIC), a-score optimization, 5 PCs, 3 discriminant functions
- **STRUCTURE**: K=3 from external software; Q-matrix integrated
- Population covariates exported for downstream GWAS

### Phase 3 — Association Mapping (`03_association_mapping.Rmd`)
- GWASpoly configured for ploidy=4
- Six models tested: additive, general, 1-dom-alt, 1-dom-ref, 2-dom-alt, 2-dom-ref
- Each trait analyzed under PCA, DAPC, and STRUCTURE corrections
- Significance assessed by Bonferroni correction (α = 0.05)
- Post-GWAS: SNP effect tables, Manhattan & QQ plots, QTL boxplots

---

## Principal Findings

| Trait | Correction | Sig. SNPs | Chromosomes |
|-------|-----------|-----------|-------------|
| Yield | PCA | 6 | 6, 7 |
| Yield | DAPC | 6 | 6, 7 |
| Yield | STRUCTURE | 6 | 6, 7 |

Six significant SNPs for yield were consistently detected on chromosomes 6 and 7, regardless of the population structure method used — supporting the robustness of these associations.

---

## Shiny Dashboard

An interactive dashboard (`R/app.R`) enables exploration of:

- Project overview & pipeline summary
- Phenotype distributions, trait correlations, top/bottom varieties
- PCA scatter plots, DAPC cluster membership, scree analysis
- Interactive Manhattan & QQ plots with model/method selectors
- Significant SNP tables and QTL boxplots

### Launch Locally
```r
install.packages(c("shiny", "shinydashboard", "plotly", "DT", "ggplot2", "dplyr", "tidyr"))
shiny::runApp("R/app.R")
```

---

## Dependencies

| Package | Role |
|---------|------|
| `GWASpoly` | Polyploid GWAS engine |
| `adegenet` | Genlight objects, PCA, DAPC |
| `lme4` | Mixed-effects modeling for BLUEs |
| `tidyverse` | Data wrangling & visualization |
| `shiny` / `shinydashboard` | Dashboard framework |
| `plotly` | Interactive plotting |
| `DT` | Interactive data tables |

---

## References

- Sharma, S.K. et al. (2018). Linkage Disequilibrium and Evaluation of Genome-Wide Association Mapping Models in Tetraploid Potato. *G3*, 8(10), 3185-3202.
- Thia, J.A. (2023). Guidelines for standardizing the application of discriminant analysis of principal components to genotype data. *Molecular Ecology Resources*, 23(3), 523-538.
- Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. *Bioinformatics*, 24(11), 1403-5.
- Bates, D. et al. (2015). Fitting Linear Mixed-Effects Models Using lme4. *Journal of Statistical Software*, 67(1), 1-48.

---

## Author

**Atharv Yeole**  
MS Bioinformatics  
2025

---

## License

MIT License — see [LICENSE](LICENSE) for details.
