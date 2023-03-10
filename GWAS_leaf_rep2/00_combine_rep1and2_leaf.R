## read data
setwd("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/")

polarity <- "negative"

if (polarity == "negative") {
    rep2 <- read.table("rep2_negative_match.csv", sep = "\t", dec = ".", header = TRUE)
    rep_feng <- read.table("all clusters_LCMS_Exp.gda", sep = "\t", dec = ".", header = TRUE)
}

## remove from rep2 the samples from seeds and QC samples
rep2 <- rep2[, !grepl(colnames(rep2), pattern = "^Sample")]
rep2 <- rep2[, !grepl(colnames(rep2), pattern = "QC")]

## rename the leaf samples to accession
ecotype_2 <- openxlsx::read.xlsx("match_accession_names_to_id.xlsx", sheet = 1) |>
    dplyr::select(sample.name, ecotype.name)
ecotype_2$sample.name <- gsub(ecotype_2$sample.name, pattern = "_pos", 
    replacement = "_Negative")
ecotype_2$sample.name <- gsub(ecotype_2$sample.name, pattern = "^X", 
    replacement = "leaf_Sample")
accession_2 <- ecotype_2$ecotype.name
names(accession_2) <- ecotype_2$sample.name
samp_cols <- colnames(rep2)[grep(colnames(rep2), pattern = "^leaf_Sample")]
colnames(rep2)[colnames(rep2) %in% samp_cols] <- as.character(accession_2[samp_cols])

## remove the columns that were not matched/not needed
rep2 <- rep2[, !is.na(colnames(rep2))]

## remove the first few lines and filter the data set
rep_feng <- rep_feng[-c(1:4), ]
rep_feng <- rep_feng[rep_feng[, "RT"] > 0.2, ]
rep_feng <- rep_feng[rep_feng[, "RT"] < 17, ]
colnames(rep_feng)[colnames(rep_feng) == "m.z"] <- "mz"

## remove the QC samples 
#rep_feng <- rep_feng[, !grepl(colnames(rep_feng), pattern = "QC")]

## rename the columns to ecotype.XYZ
ecotype_feng <- openxlsx::read.xlsx("data information.xlsx", sheet = 1) |>
    dplyr::select(sample.name, accession.name) 
ecotype_feng$sample.name <- paste0("X", ecotype_feng$sample.name)
accession_feng <- ecotype_feng$accession.name
names(accession_feng) <- ecotype_feng$sample.name
samp_cols <- colnames(rep_feng)[grep(colnames(rep_feng), pattern = "^X")]
colnames(rep_feng)[colnames(rep_feng) %in% samp_cols] <- as.character(accession_feng[samp_cols])

## remove the columns that were not matched/not needed
rep_feng <- rep_feng[, !is.na(colnames(rep_feng))]
#rep_feng <- rep_feng[, !grepl(colnames(rep_feng), pattern = "extraction QC")]
#rep_feng <- rep_feng[, !grepl(colnames(rep_feng), pattern = "running QC")]
rep_feng <- rep_feng[, !grepl(colnames(rep_feng), pattern = "ecotype[.]6909[.]")]

## get the internal standard: RT 6.96, m/z 431.097847
int_std <- rep_feng[rep_feng$mz > 431.08 & rep_feng$mz < 431.11 & rep_feng$RT >6.8 & rep_feng$RT < 7.0 , ]
saveRDS(int_std, file = "internal_standard_feng.RDS")

rt.dev.range <- c(-0.8, 0.8)
mz.dev.range <- c(-0.01, 0.01)

res <- vector("list", nrow(rep2))
for (i in 1:nrow(rep2)) {
    mz.diff <- rep2$mz[i] - rep_feng$mz
    rt.diff <- rep2$RT[i] - rep_feng$RT
    j <- which(rt.diff > rt.dev.range[1] & rt.diff < rt.dev.range[2] &
        mz.diff > mz.dev.range[1] & mz.diff < mz.dev.range[2])
    res[[i]] <- if (length(j) == 0) NA else j
}

## obtain values only
if (polarity == "negative") {
    rep2_values <- rep2[, which(colnames(rep2) == "ecotype.7283"):which(colnames(rep2) == "ecotype.6188")]
    rep_feng_values <- rep_feng[, which(colnames(rep_feng) == "ecotype.173"):which(colnames(rep_feng) == "ecotype.6909")]
    .cols <- intersect(colnames(rep2_values), colnames(rep_feng_values))
    rep2_values <- rep2_values[, .cols]
    rep_feng_values <- rep_feng_values[, .cols]
}

## convert to matrix
rep2_values <- as.matrix(rep2_values)
mode(rep2_values) <- "numeric"
rep2_values[rep2_values == 0] <- NA
rep_feng_values <- as.matrix(rep_feng_values)
mode(rep_feng_values) <- "numeric"

## check that we have the same column names
all(colnames(rep2_values) == colnames(rep_feng_values))

## filter: only take the features that are present in both replicates, for features in rep1 that 
## were mapped multiple times 
## go through res (rows in rep1) and keep the corresponding row in rep2 that has the highest covariance term
## rep1
cov_res <- rep(NaN, length(res))
res_cor <- vector("list", length(res))
for (i in 1:length(res)) {
    inds <- res[[i]]
    if (!all(is.na(inds))) {
        if (length(inds) > 1) {
            covs <- lapply(inds, function(inds_i) {
                tryCatch(
                    cov(log2(rep2_values[i, ] + 1), log2(rep_feng_values[inds_i, ] + 1), use = "complete.obs"),
                    error = function(e) NA)
            }) |> unlist()
            cov_max <- covs[which.max(covs)]
            if (length(cov_max) > 0) {
                cov_res[i] <- cov_max
                res_cor[[i]] <- res[[i]][which.max(covs)]      
            } else {
                cov_res[i] <- NaN
                res_cor[[i]] <- NaN
            }
        } else {
            covs <- tryCatch(expr = {
                cov(log2(rep2_values[i, ] + 1), log2(rep_feng_values[inds,] + 1), use = "complete.obs")
            }, error = function(e) NA)
            if (is.na(covs)) {
                cov_res[i] <- NaN
                res_cor[[i]] <- NaN
            } else {
                cov_res[i] <- covs
                res_cor[[i]] <- res[[i]]
            }
        }
    } else {
        cov_res[i] <- NaN
        res_cor[[i]] <- NaN
    }
    
}

res_cor <- unlist(res_cor)

## set features whose features is equal or smaller than 0 to NaN
res_cor[cov_res <= 0] <- NaN
cov_res[cov_res <= 0] <- NaN

rep2_cor <- rep2[!is.na(res_cor), ]
rep_feng_cor <- rep_feng[res_cor[!is.na(res_cor)], ]

## add the relation of the mapping
rep2_cor <- cbind(rep2_cor, 
    mapping_rep2_rep_feng = paste(rep2_cor[, "Name"], rep_feng[res_cor[!is.na(res_cor)], "Name"], sep = "/"))
rep_feng_cor <- cbind(rep_feng_cor, 
    mapping_rep2_rep_feng = paste(rep2[!is.na(res_cor), "Name"], rep_feng_cor[, "Name"], sep = "/"))

## add the covariance with the mapped feature
rep2_cor <- cbind(rep2_cor, cov = cov_res[!is.na(cov_res)])
rep_feng_cor <- cbind(rep_feng_cor, cov = cov_res[!is.na(cov_res)])

## add RT deviance 
rep2_cor <- cbind(rep2_cor, rt_dev = rep2_cor[, "RT"] - rep_feng_cor[, "RT"])
rep_feng_cor <- cbind(rep_feng_cor, rt_dev = rep2_cor[, "RT"] - rep_feng_cor[, "RT"])

## add m/z deviance
rep2_cor <- cbind(rep2_cor, mz_dev = rep2_cor[, "mz"] - rep_feng_cor[, "mz"])
rep_feng_cor <- cbind(rep_feng_cor, mz_dev = rep2_cor[, "mz"] - rep_feng_cor[, "mz"])

## save
##write.table(rep2_cor, file = paste0("./rep2_", polarity, "_match_with_feng.csv"), sep = "\t", dec = ".", quote = FALSE)
write.table(rep_feng_cor, file = paste0("./rep_feng_", polarity, "_match.csv"), sep = "\t", dec = ".", quote = FALSE)

