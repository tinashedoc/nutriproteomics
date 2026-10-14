# Nutriproteomics GWAS

Genome-wide association analysis of four nutrient pattern proteomic scores
(plant-driven, animal-driven, vitamin C/sugar-driven and retinol/vitamin B12-driven)
using fastGWA-MLM in [GCTA](https://yanglab.westlake.edu.cn/software/gcta/).

## Pipeline

| Step | Script | Output |
|------|--------|--------|
| 1. Transform phenotypes | `scripts/01_transform_phenotypes.R` | `pheno/pheno{plant,animal,vitc,retinol}.txt`, `pheno/qcovar.txt`, `pheno/covar.txt` |
| 2. Run GWAS | `scripts/02_run_fastgwa.sh` | `gwas/sparse_grm.*`, `gwas/<pheno>.fastGWA` |
| 3. Summarise | `scripts/03_gwas_summary.R` | genomic inflation (lambda), `gwas/<pheno>_qq.png`, `gwas/<pheno>_manhattan.png` |

Requirements: R (base packages only) and GCTA 1.94.1.

```bash
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
unzip gcta-1.94.1-linux-kernel-3-x86_64.zip
```

## 1. Transforming the proteomic scores

Each score is regressed on age, sex and the first 10 genetic principal components.
The residuals are inverse rank normalised and written as GCTA phenotype files
(`FID IID value`, no header).

```bash
Rscript scripts/01_transform_phenotypes.R DATASET.tsv pheno
```

`DATASET.tsv` (tab- or whitespace-separated; `.csv` also accepted) needs one row per
individual with these columns:

| Column | Description |
|--------|-------------|
| `FID`, `IID` | Family and individual IDs, matching the PLINK `.fam` file |
| `age` | Age in years |
| `sex` | Sex (treated as categorical) |
| `PC1` … `PC10` | Genetic principal components |
| `plantdriven`, `animaldriven`, `vitaminCsugar`, `retinolvitB12` | Nutrient pattern proteomic scores |

Column names, the number of PCs and the covariate set are defined at the top of the script.
Missing values are allowed: `lm(..., na.action = na.exclude)` keeps the residuals
aligned with their individuals, and those individuals get `NA` in the phenotype file.

The core of the transformation, for one score:

```R
fit <- lm(plantdriven ~ age + sex + PC1 + ... + PC10, data = DATASET, na.action = na.exclude)
DATASET$residualsplant <- residuals(fit)
DATASET$phenoplant <- qnorm((rank(DATASET$residualsplant, na.last = "keep") - 0.5) /
                            sum(!is.na(DATASET$residualsplant)))
```

## 2. GWAS with fastGWA-MLM

The genotype data must be in PLINK binary format. Configure the run with environment variables:

```bash
GCTA=/path/to/gcta-1.94.1 \
BFILE=/path/to/all_imputed_map_qc \
THREADS=20 \
GRM_SNPS=pruned_genotyped.snplist \
bash scripts/02_run_fastgwa.sh
```

| Variable | Default | Meaning |
|----------|---------|---------|
| `GCTA` | `gcta-1.94.1` | GCTA executable |
| `BFILE` | required | PLINK genotype prefix |
| `PHENO_DIR` | `pheno` | Output of step 1 |
| `OUT_DIR` | `gwas` | Results directory |
| `THREADS` | `20` | Threads |
| `GRM_SNPS` | none | Optional SNP list for the GRM; LD-pruned, genotyped SNPs are recommended over all imputed SNPs |
| `GRM_MAF` | `0.01` | MAF filter for SNPs in the GRM |
| `SPARSE_CUTOFF` | `0.05` | Relatedness cutoff for the sparse GRM |
| `GWAS_MAF` | `0.01` | MAF filter for tested SNPs |
| `GWAS_GENO` | `0.05` | Maximum per-SNP missingness for tested SNPs |

The script builds a sparse GRM once and then runs fastGWA-MLM for each phenotype.
Age, PCs (`--qcovar`) and sex (`--covar`) are included again at the GWAS stage.
Inverse normalising residuals and then testing without covariates can lose power and
bias effect estimates (Sofer et al., 2019, *Genet Epidemiol* 43:263–275).

Filter imputed variants on imputation quality (e.g. INFO/R² ≥ 0.3) before creating the
PLINK files, or use dosages directly with GCTA's `--bgen`/`--mbgen` options.

## 3. Summarising results

```bash
Rscript scripts/03_gwas_summary.R gwas
```

This reports the genomic inflation factor (lambda) and the number of genome-wide significant
SNPs (P < 5e-8) for each phenotype, and saves QQ and Manhattan plots.
