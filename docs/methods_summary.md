# Analytical Methods

## Experimental Design
- **Organism**: Autotetraploid potato (*Solanum tuberosum*), ploidy = 4
- **Diversity panel**: 282 varieties
- **Marker platform**: Illumina Infinium 8K SNP array
- **Field trials**: 2016-2018 at Nafferton Farm, Newcastle University
- **Treatments**: Organic and Conventional systems, 2 replicates per system

## Phenotypic Traits
| Trait | Unit | Seasons |
|-------|------|---------|
| Yield | kg | 2016-2018 |
| Max Max Height | m | 2017-2018 |
| Max Average Height | m | 2017-2018 |
| Max Canopy Area | m² | 2017-2018 |
| Max Canopy Volume | m³ | 2017-2018 |

## Workflow

### Genotype Quality Control
- Markers with ≥20% missing data excluded (per Sharma et al., 2018)
- PIC > 0.4 threshold applied for STRUCTURE/DAPC comparability

### BLUEs Estimation
- Linear mixed model: `Trait ~ Genotype + (1|System) + (1|Year)`
- REML estimation via lme4
- Model comparison using log-likelihood ratio tests, AIC, BIC
- Missing phenotypes imputed from same-variety replicates

### Population Structure
- **PCA**: glPca() on genlight object; 5 PCs retained based on scree analysis
- **DAPC**: K-means (K=4 from BIC), a-score optimization, 5 PCs, 3 discriminants
- **STRUCTURE**: K=3 via external software; Q-matrix imported as covariates

### GWAS (GWASpoly)
- Genetic models: additive, general, 1-dom-alt, 1-dom-ref, 2-dom-alt, 2-dom-ref
- Population corrections: PCA, DAPC, STRUCTURE (run independently)
- Kinship matrix computed internally
- Multiple testing: Bonferroni correction at α = 0.05

## Results Summary
- 6 significant SNPs detected for Yield on chromosomes 6 and 7
- Consistent across all population structure correction methods
- General model yielded strongest associations under DAPC correction
