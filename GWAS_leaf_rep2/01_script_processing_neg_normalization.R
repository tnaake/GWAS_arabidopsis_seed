################################################################################
setwd("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/")

## load libraries
library(SummarizedExperiment)
library(MatrixQCvis)

## read the data
met <- read.table("rep_feng_negative_match.csv", sep = "\t", quote = "")
meta_data <- openxlsx::read.xlsx("running days weight information.xlsx")
int_std <- readRDS("internal_standard_feng.RDS")

## create SummarizedExperiment object

## assay 
.cols <- which(colnames(met) == "ecotype.173"):which(colnames(met) == "ecotype.6909")
a <- met[, .cols] |>
    as.matrix()
mode(a) <- "numeric"
a <- a[, !grepl(colnames(a), pattern = "QC")]

## rowData
rD <- met[, c("Name", "Mass", "mz", "RT", "mapping_rep2_rep_feng", "cov", "rt_dev", "mz_dev")]

## colData
cD <- meta_data[meta_data[["accession_name"]] %in% colnames(a), ]
cD <- cD[match(colnames(a), cD[["accession_name"]]), ]
rownames(cD) <- cD[["accession_name"]]

## create SummarizedExperiment
se <- SummarizedExperiment(assays = a, rowData = rD, colData = cD)
## check by shinyQC(se)

## quality check by MatrixQCvis, remove 7310 and 7062
se <- MatrixQCvis:::selectSampleSE(se, c("ecotype.7310", "ecotype.7062"), 
    mode = "exclude")
se <- se[, !is.na(se$weight_g)]
## check by shinyQC(se)

## remove rows that contain less than five measured values
se <- se[apply(assay(se), 1, function(rows_i) sum(!is.na(rows_i))) >= 10, ]
## check by shinyQC(se)

## 1) divide by the weight
se <- sweep(x = assay(se), MARGIN = 2, FUN = "/", STATS = se$weight_g) |>
    MatrixQCvis:::updateSE(se = se)

## 2) perform log2 transformation
se <- assay(se) |>
    transformAssay(method = "log2") |>
    MatrixQCvis:::updateSE(se = se)
## check by shinyQC(se)

## 3) perform batch correction
se <- se |>
    batchCorrectionAssay(method = "removeBatchEffect (limma)", 
        batchColumn = "batch") |>
    MatrixQCvis:::updateSE(se = se)
## check by shinyQC(se)

## looking at dimension reduction plots, there seems a batch effect of unknown
## origin, determine the samples that belong to the two clusters (stored in 
## "unknown_factor"), remove some samples where is is unsure to which cluster
## they belong to
se <- MatrixQCvis:::selectSampleSE(se, 
    c("ecotype.6962", "ecotype.6911", "ecotype.8376", "ecotype.6940",
        "ecotype.8271", "ecotype.7098", "ecotype.6933", "ecotype.242"), 
    mode = "exclude")
se_batch <- se |>
    batchCorrectionAssay(method = "removeBatchEffect (limma)",
        batchColumn = "unknown_cofactor") |>
    MatrixQCvis:::updateSE(se = se)
## check by shinyQC(se)
## check by shinyQC(se_batch)

rownames(rowData(se)) <- rowData(se)$Name
rownames(se) <- rowData(se)$Name
rownames(rowData(se_batch)) <- rowData(se_batch)$Name
rownames(se_batch) <- rowData(se_batch)$Name

saveRDS(se, file="normalization_leaf_feng.RDS")
saveRDS(se_batch, file="normalization_leaf_feng_batch_corrected.RDS")
