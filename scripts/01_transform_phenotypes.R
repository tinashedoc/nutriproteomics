#!/usr/bin/env Rscript
# Transform the four nutrient pattern proteomic scores into GWAS phenotypes.
#
# For each score: regress on age, sex and genetic PCs, inverse rank normalise
# the residuals, and write a GCTA phenotype file (FID IID value, no header).
# Also writes the covariate files used again at the GWAS stage.
#
# Usage: Rscript 01_transform_phenotypes.R <input.tsv|csv> [output_dir]
#
# The input needs one row per individual with columns FID and IID (matching
# the PLINK .fam file), age, sex, PC1..PC<N_PCS> and the score columns below.

## ---- Configuration -------------------------------------------------------

# Phenotype label (used in output file names) -> score column in the input.
SCORES <- c(
  plant   = "plantdriven",
  animal  = "animaldriven",
  vitc    = "vitaminCsugar",
  retinol = "retinolvitB12"
)

# Number of genetic principal components to adjust for.
N_PCS <- 10

QCOVARS <- c("age", paste0("PC", seq_len(N_PCS)))  # quantitative
COVARS  <- c("sex")                                # categorical

## ---- Helpers ---------------------------------------------------------------

# Rank-based inverse normal transformation; NAs are kept in place.
inverse_normal <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

write_gcta <- function(df, path) {
  write.table(df, path, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE, na = "NA")
}

## ---- Main --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 01_transform_phenotypes.R <input.tsv|csv> [output_dir]")
}
input   <- args[1]
out_dir <- if (length(args) >= 2) args[2] else "pheno"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sep <- if (grepl("\\.csv$", input, ignore.case = TRUE)) "," else ""
DATASET <- read.table(input, header = TRUE, sep = sep,
                      stringsAsFactors = FALSE, check.names = FALSE)

required <- c("FID", "IID", QCOVARS, COVARS, SCORES)
missing  <- setdiff(required, names(DATASET))
if (length(missing) > 0) {
  stop("Input is missing columns: ", paste(missing, collapse = ", "))
}

DATASET$sex <- factor(DATASET$sex)

for (label in names(SCORES)) {
  score   <- SCORES[[label]]
  formula <- reformulate(c(QCOVARS, COVARS), response = score)

  # na.exclude pads residuals with NA so they stay aligned with DATASET rows.
  fit <- lm(formula, data = DATASET, na.action = na.exclude)

  resid_col <- paste0("residuals", label)
  pheno_col <- paste0("pheno", label)
  DATASET[[resid_col]] <- residuals(fit)
  DATASET[[pheno_col]] <- inverse_normal(DATASET[[resid_col]])

  out <- file.path(out_dir, paste0(pheno_col, ".txt"))
  write_gcta(DATASET[, c("FID", "IID", pheno_col)], out)
  message(sprintf("%-8s %s ~ %d covariates, n = %d -> %s",
                  label, score, length(c(QCOVARS, COVARS)),
                  sum(!is.na(DATASET[[pheno_col]])), out))
}

write_gcta(DATASET[, c("FID", "IID", QCOVARS)], file.path(out_dir, "qcovar.txt"))
write_gcta(DATASET[, c("FID", "IID", COVARS)],  file.path(out_dir, "covar.txt"))
message("Covariate files written to ", out_dir)
