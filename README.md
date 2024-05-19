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
DATASET$phenoanimal <- qnorm((rank(DATASET$residuals,na.last="keep") - 0.5)/sum(!is.na(DATASET$residuals)))
DATASET$phenovitC <- qnorm((rank(DATASET$residuals,na.last="keep") -0.5)/sum(!is.na(DATASET$residuals)))
DATASET$phenovitC <- qnorm((rank(DATASET$residuals,na.last="keep") -0.5)/sum(!is.na(DATASET$residuals)))
