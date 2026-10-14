#!/usr/bin/env bash
# Run fastGWA-MLM (GCTA) for the four nutrient pattern proteomic phenotypes.
#
# Run 01_transform_phenotypes.R first. Settings can be overridden with
# environment variables, e.g.
#   GCTA=/path/to/gcta-1.94.1 BFILE=/data/all_imputed_map_qc ./02_run_fastgwa.sh
set -euo pipefail

GCTA=${GCTA:-gcta-1.94.1}           # GCTA executable
BFILE=${BFILE:?Set BFILE to the PLINK genotype prefix (without .bed/.bim/.fam)}
PHENO_DIR=${PHENO_DIR:-pheno}       # output directory of 01_transform_phenotypes.R
OUT_DIR=${OUT_DIR:-gwas}
THREADS=${THREADS:-20}
GRM_SNPS=${GRM_SNPS:-}              # optional SNP list (e.g. LD-pruned, genotyped) for the GRM
GRM_MAF=${GRM_MAF:-0.01}
SPARSE_CUTOFF=${SPARSE_CUTOFF:-0.05}
GWAS_MAF=${GWAS_MAF:-0.01}
GWAS_GENO=${GWAS_GENO:-0.05}        # max per-SNP missingness

PHENOS=(plant animal vitc retinol)

mkdir -p "$OUT_DIR"
GRM="$OUT_DIR/sparse_grm"

# Sparse genetic relationship matrix from common (optionally pruned) SNPs.
if [[ ! -f "$GRM.grm.id" ]]; then
  grm_args=(--bfile "$BFILE" --maf "$GRM_MAF" --make-grm
            --sparse-cutoff "$SPARSE_CUTOFF" --thread-num "$THREADS" --out "$GRM")
  [[ -n "$GRM_SNPS" ]] && grm_args+=(--extract "$GRM_SNPS")
  "$GCTA" "${grm_args[@]}"
fi

# Covariates are included again at the GWAS stage (fully adjusted two-stage).
for pheno in "${PHENOS[@]}"; do
  "$GCTA" --bfile "$BFILE" \
    --grm-sparse "$GRM" \
    --fastGWA-mlm \
    --pheno "$PHENO_DIR/pheno$pheno.txt" \
    --qcovar "$PHENO_DIR/qcovar.txt" \
    --covar "$PHENO_DIR/covar.txt" \
    --maf "$GWAS_MAF" \
    --geno "$GWAS_GENO" \
    --thread-num "$THREADS" \
    --out "$OUT_DIR/$pheno"
done
