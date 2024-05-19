# Nutriproteomics

**Transformation of the proteomic markers**

Regress nutrient pattern proteomic scores with appropriate covariates, including age, sex, and genetic principal components. Inverse rank normalize the residuals and uses the transformed phenotype as the outcome variable or phenotype for  GWAS. The following R commands can be used to transform your trait.Regress nutrient pattern proteomic scores with appropriate covariates, including age, sex, and genetic principal components. Inverse rank normalize the residuals and uses the transformed phenotype as the outcome variable or phenotype for  GWAS. The following R commands can be used to transform your trait.

```R
#Regress nutrient proteomic scores with select covariates (separately for the four nutrient pattern proteomic score)

plantdrivenlm <-lm(plantdriven ~ PC1 + PC2 + PC3 + PC4 + PC5 + age, data=DATASET)
animaldrivenlm <-lm(animaldriven ~ PC1 + PC2 + PC3 + PC4 + PC5 + age, data=DATASET)
vitamicsugardrivenlm <- lm(vitaminCsugar ~ PC1 + PC2 + PC3 + PC4 + PC5 + age +sex,data=DATASET)
retinolvitB12drivenlm <- lm(vitaminCsugar ~ PC1 + PC2 + PC3 + PC4 + PC5 + age +sex,data=DATASET)

#Extract residuals and add them to your dataset

DATASET$residualsplant = plantdrivenlm$residuals
DATASET$residualsanimal = animaldrivenlm$residuals
DATASET$residualsvitC = vitamicsugardrivenlm$residuals
DATASET$residualsretinol = retinolvitB12drivenlm$residuals

#Inverse rank normalizes the residuals and renames them as pheno for GWAS analysis

DATASET$phenoplant <- qnorm((rank(DATASET$residualsplant,na.last="keep") -0.5)/sum(!is.na(DATASET$residualsplant)))
DATASET$phenoanimal <- qnorm((rank(DATASET$residualsanimal,na.last="keep") - 0.5)/sum(!is.na(DATASET$residualsanimal)))
DATASET$phenovitC <- qnorm((rank(DATASET$residualsvitc,na.last="keep") -0.5)/sum(!is.na(DATASET$residualsvitc)))
DATASET$phenoretinol <- qnorm((rank(DATASET$residualsretinol,na.last="keep") -0.5)/sum(!is.na(DATASET$residualsretinol)))
```

**GWAS analysis using FastGWA**

```Bash

#Download the Fastgwas executable from link below

wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
unzip gcta-1.94.1-linux-kernel-3-x86_64.zip

#Make sure you genotype data is in plink formart and use the Fastgwas to run the GWAS as indicated below

#Define the location of the Fastgwas executable as SOFT and the genotype data as DATA
DATA=/dataE/AWIGenGWAS/shared/imputed_data_plink

SOFT=/home/tinashe/GCTA/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1

#Make a sparse genetic related matrix

$SOFT --bfile $DATA/all_imputed_map_qc --make-grm --sparse-cutoff 0.05 --thread-num 20 --out assp_grm

#Run the GWAS as shown below

$SOFT --bfile $DATA/all_imputed_map_qc --grm-sparse assp_grm --fastGWA-mlm --pheno phenoplant --thread-num 20 --out plant

$SOFT --bfile $DATA/all_imputed_map_qc --grm-sparse assp_grm --fastGWA-mlm --pheno phenoanimal --thread-num 20 --out animal

$SOFT --bfile $DATA/all_imputed_map_qc --grm-sparse assp_grm --fastGWA-mlm --pheno phenovitc --thread-num 20 --out vitc

$SOFT --bfile $DATA/all_imputed_map_qc --grm-sparse assp_grm --fastGWA-mlm --pheno phenoretinol --thread-num 20 --out retinol











