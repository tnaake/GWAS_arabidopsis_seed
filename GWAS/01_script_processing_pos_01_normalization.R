library(gdps)
memory.limit(56000)
###############################################################################
setwd("/home/naake/01_GWAS/")
met1 <- read.table("rep1_positive_match.csv", sep = "\t", quote = "")
met2 <- read.table("rep2_positive_match.csv", sep = "\t", quote = "")

## match sample names and table with weights
## load excel sheet
library(xlsx)
weights_seeds1 <- read.xlsx("GWAS_accession_I_samplelist_weights_bioRep1.xlsx", sheetIndex = 1, startRow = 2)
weights_seeds2 <- read.xlsx("GWAS_accession_I_samplelist_weights_bioRep2.xlsx", sheetIndex = 1, startRow = 2)
colnames(weights_seeds1)[grep(colnames(weights_seeds1), pattern = "weight")] <- "weight."
weights_seeds1[, "weight."] <- weights_seeds1[, "weight."] / 1000 ## in mg
weights_seeds2[, "weight."] <- weights_seeds2[, "weight."] / 1000 ## in mg
weights_leaf <- read.xlsx("GWAS_leaf_Si_control_weight.xlsx", sheetIndex = 1, startRow = 1) # in mg

divide_by_weight <- function(met, weight, tissue = "seed") {
    
    met <- as.matrix(met)
    mode(met) <- "numeric"
    
    ## rename column with weigths
    colnames(weight)[grep(colnames(weight), pattern = "weight")] <- "weight"
    
    ind_cut <- unlist(lapply(strsplit(colnames(met), split = "Sample"), "[", 2))    
    ind_cut <- unlist(lapply(strsplit(ind_cut, split = "_"), "[", 1))
    
    ## determine the row indices in weights_seeds to take for finding the 
    ## corresponding weight
    if (tissue=="seed") ind_to_sweep <- ifelse(ind_cut %in% 
        weight[, "Internal.ID"], match(ind_cut, weight[, "Internal.ID"]), FALSE)
    if (tissue=="leaf") ind_to_sweep <- ifelse(ind_cut %in% 
        weight[, "NO"], match(ind_cut, weight[, "NO"]), FALSE)
    ## create vector that contains the weights 
    ## (and matches in order to the columns in met_new)
    weights_dup <- weight[ind_to_sweep, "weight"]

    ## divide by weights
    met <- sweep(x = met, MARGIN = 2, STATS = weights_dup, FUN = "/")
    return(met)   
}


cols_seed1_pos <- which(colnames(met1) == "Sample1_Positive"):which(colnames(met1) == "Sample338_Positive")
cols_leaf_QC_pos <- which(colnames(met2) == "leaf_QC1__begin_Positive"):which(colnames(met2) == "leaf_QC6_end_Positive")
cols_leaf_pos <- which(colnames(met2) == "leaf_Sample1_Positive"):which(colnames(met2) == "leaf_Sample336_Positive")
cols_seed2_QC_pos <- which(colnames(met2) == "QC1_Positive"):which(colnames(met2) == "QC6_end_Positive")
cols_seed2_pos <- which(colnames(met2) == "Sample1_Positive"):which(colnames(met2) == "Sample338_Positive")

weight_seed2_QC <- mean(c(0.466, 0.481, 0.483, 0.520, 0.556, 0.608))

## samples that were used to create mix
weight_leaf_QC <- c(10, 11, 23, 29, 50,  51, 53, 56, 58, 59, 65, 66, 68, 69, 
    71, 72,  78, 80, 81, 96, 101, 105, 106, 109, 110, 114, 117, 124, 125, 126, 
    127, 128, 129, 130, 132, 135, 137, 138, 139, 144, 145, 146,  147, 148, 150, 
    152, 154, 155, 156, 157, 158, 159, 160, 164, 166, 170, 176, 177, 178, 181, 
    183, 184, 187, 193, 196, 200, 203, 204, 210, 213, 216, 225, 233, 253, 255, 
    257, 258,  270, 260, 272, 283, 284, 293, 302, 304, 305, 309, 311, 313, 314, 
    318, 321, 326, 327, 328, 329, 330, 332, 335) 
weight_leaf_QC <- mean(weights_leaf[weights_leaf[,"NO"] %in% weight_leaf_QC, "weight"]) 

met1_dw <- met1
met2_dw <- met2

## divide by weight
met1_dw[, cols_seed1_pos] <- divide_by_weight(met1_dw[, cols_seed1_pos], weight = weights_seeds1)
met2_dw[, cols_leaf_QC_pos] <- met2_dw[, cols_leaf_QC_pos] / weight_leaf_QC
met2_dw[, cols_leaf_pos] <- divide_by_weight(met2_dw[, cols_leaf_pos], weight = weights_leaf, tissue = "leaf")
met2_dw[, cols_seed2_QC_pos] <- met2_dw[, cols_seed2_QC_pos] / weight_seed2_QC
met2_dw[, cols_seed2_pos] <- divide_by_weight(met2_dw[, cols_seed2_pos], weight = weights_seeds2)

## remove double features from met2
met2_dw_rem <- met2_dw[!duplicated(met2_dw[, "Name"]), ]

boxplot(sqrt(log2(met1_dw[, 2:314] + 1)), use.cols = TRUE)
boxplot(sqrt(log2(met2_dw_rem[, 2:309] + 1)), use.cols = TRUE) ## remove leaf_QC6__begin_positive
boxplot(log2(met2_dw_rem[, 310:631] + 1), use.cols = TRUE)

## remove 124, 130, 174, 237 from seed1
met1_dw <- met1_dw[, -which(colnames(met1_dw) %in% 
    c("Sample124_Positive", "Sample130_Positive", "Sample174_Positive", "Sample237_Positive"))]

## normalize
normalize <- function(met, batch_list, qc = NULL) {
    
    batch1 <- which(colnames(met) == batch_list[[1]][1]):which(colnames(met) == batch_list[[1]][2])
    if (!is.null(qc)) batch1 <- c(which(colnames(met) %in% qc[[1]]), batch1)
    batch1_l <- length(batch1)
    batch2 <- which(colnames(met) == batch_list[[2]][1]):which(colnames(met) == batch_list[[2]][2])
    if (!is.null(qc)) batch2 <- c(which(colnames(met) %in% qc[[2]]), batch2)
    batch2_l <- length(batch2)
    batch3 <- which(colnames(met) == batch_list[[3]][1]):which(colnames(met) == batch_list[[3]][2])
    if (!is.null(qc)) batch3 <- c(which(colnames(met) %in% qc[[3]]), batch3)
    batch3_l <- length(batch3)
    batch4 <- which(colnames(met) == batch_list[[4]][1]):which(colnames(met) == batch_list[[4]][2])
    if (!is.null(qc)) batch4 <- c(which(colnames(met) %in% qc[[4]]), batch4)
    batch4_l <- length(batch4)
    batch5 <- which(colnames(met) == batch_list[[5]][1]):which(colnames(met) == batch_list[[5]][2])
    if (!is.null(qc)) batch5 <- c(which(colnames(met) %in% qc[[5]]), batch5)
    batch5_l <- length(batch5)
    batch6 <- which(colnames(met) == batch_list[[6]][1]):which(colnames(met) == batch_list[[6]][2])
    if (!is.null(qc)) batch6 <- c(which(colnames(met) %in% qc[[6]]), batch6)
    batch6_l <- length(batch6)
    
    col_day <- c(rep("blue", batch1_l), rep("red", batch2_l), 
        rep("yellow", batch3_l), rep("green", batch4_l), 
        rep("grey", batch5_l), rep("cyan", batch6_l))
    day <- c(rep(1, batch1_l), rep(2, batch2_l), rep(3, batch3_l), 
             rep(4, batch4_l), rep(5, batch5_l), rep(6, batch6_l))
    
    inds_batch <- c(batch1, batch2, batch3, batch4, batch5, batch6)
    smp <- data.frame(sample = colnames(met[, inds_batch]), 
        colour = col_day, day = day)
    rownames(smp) <- smp[, 1]
    
    ## log2 normalization
    mat_log <- log2(met[, inds_batch] + 1)
    
    ## Day normalization (DN)
    uniqueMD <- unique(smp[, 3])
    mat_log_DNtemp <- mat_log

    for (MD in uniqueMD) {
        filtMD <- smp[, 3] == MD
        subdata <- mat_log[, filtMD]
        metab.med <- apply(subdata, 1, median, na.rm = TRUE) 
        mat_log_DNtemp[, filtMD] <- sweep(subdata, 1, metab.med, FUN = "-")
    }
    
    ## add median per metabolite of whole experiment.
    metab.med.all <- apply(mat_log, 1, median, na.rm = TRUE) 
    mat_log_DN <- mat_log_DNtemp + metab.med.all

    mat_log_DMN <- limma::removeBatchEffect(x = mat_log, batch = smp[, 3])
    
    rownames(smp) <- smp[, 1]
    mat_log <- as.matrix(mat_log)
    mat_log_DN <- as.matrix(mat_log_DN)
    mat_log_DMN <- as.matrix(mat_log_DMN)
    
    mode(mat_log) <- "numeric"
    mode(mat_log_DN) <- "numeric"
    mode(mat_log_DMN) <- "numeric"
    
    final <- list(smp = smp, mat_log = mat_log, mat_log_DN = mat_log_DN, 
        mat_log_DMN = mat_log_DMN)

    return(final)
}

## met1 
batch_list_met1 <- list(c("Sample1_Positive", "Sample55_Positive"), 
    c("Sample56_Positive", "Sample112_Positive"), 
    c("Sample113_Positive", "Sample168_Positive"), 
    c("Sample169_Positive", "Sample225_Positive2"), 
    c("Sample226_Positive", "Sample283_Positive"), 
    c("Sample284_Positive", "Sample338_Positive"))

## met2 seed
batch_list_seed_met2 <- list(c("Sample1_Positive", "Sample55_Positive"), 
    c("Sample56_Positive", "Sample112_Positive"), 
    c("Sample113_Positive", "Sample168_Positive"), 
    c("Sample169_Positive", "Sample225_Positive"), 
    c("Sample226_Positive", "Sample283_Positive"), 
    c("Sample284_Positive", "Sample336_Positive"))
qc_seed_met2 <- list(c("QC1_Positive"), c("QC2_end_Positive"), 
    c("QC3_begin_Positive", "QC3_end_Positive"), 
    c("QC4_begin_Positive", "QC4_end_Positive"), 
    c("QC5_begin_Positive", "QC5_end_Positive"), 
    c("QC6_begin_Positive", "QC6_end_Positive"))

## met2 leaf
batch_list_leaf_met2 <- list(c("leaf_Sample1_Positive", "leaf_Sample57_Positive"), 
    c("leaf_Sample58_Positive", "leaf_Sample118_Positive"), 
    c("leaf_Sample120_Positive", "leaf_Sample166_Positive"),
    c("leaf_Sample170_Positive", "leaf_Sample226_Positive"), 
    c("leaf_Sample227_Positive", "leaf_Sample285_Positive"), 
    c("leaf_Sample286_Positive", "leaf_Sample336_Positive"))
qc_leaf_met2 <- list(c("leaf_QC1__begin_Positive", "leaf_QC1__end_Positive"), 
    c("leaf_QC2__begin_Positive", "leaf_QC2__end_Positive"), 
    c("leaf_QC3__begin_Positive", "leaf_QC3__end_Positive"), 
    c("leaf_QC4__begin_Positive", "leaf_QC4__end_Positive"), 
    c("leaf_QC5__begin_Positive", "leaf_QC5__end_Positive"), 
    c("leaf_QC6__begin_Positive", "leaf_QC6__end_Positive"))

met1_norm <- normalize(met1_dw, batch_list = batch_list_met1)
met2_seed_norm <- normalize(met2_dw_rem, batch_list = batch_list_seed_met2, qc = qc_seed_met2)
met2_leaf_norm <- normalize(met2_dw_rem, batch_list = batch_list_leaf_met2, qc = qc_leaf_met2)

save(met1_norm, met2_seed_norm, met2_leaf_norm, file = "normalization_seed_leaf_dataset_two_replicates.RData")


## PCA
library(pcaMethods)
mat_log_seed1 <- as.matrix(met1_norm$mat_log)
mat_log_DN_seed1 <- as.matrix(met1_norm$mat_log_DN)
mat_log_DMN_seed1 <- as.matrix(met1_norm$mat_log_DMN)

mat_log_seed2 <- as.matrix(met2_seed_norm$mat_log)
mat_log_DN_seed2 <- as.matrix(met2_seed_norm$mat_log_DN)
mat_log_DMN_seed2 <- as.matrix(met2_seed_norm$mat_log_DMN)

mat_log_leaf2 <- as.matrix(met2_leaf_norm$mat_log)
mat_log_DN_leaf2 <- as.matrix(met2_leaf_norm$mat_log_DN)
mat_log_DMN_leaf2 <- as.matrix(met2_leaf_norm$mat_log_DMN)

rpcaNN_seed1 <- pca(t(mat_log_seed1), method = "ppca", 
    center = TRUE, scale = 'none', nPcs = 5, seed = 3455)
rpcaDN_seed1 <- pca(t(mat_log_DN_seed1), method = "ppca", 
    center = TRUE, scale = 'none', nPcs = 5, seed = 3455, maxIterations = 1500)
rpcaDMN_seed1 <- pca(t(mat_log_DMN_seed1), method = "ppca", 
    center = TRUE, scale = 'none', nPcs = 5, seed = 3455, maxIterations = 1500)

rpcaNN_seed2 <- pca(t(mat_log_seed2), method = "ppca", 
    center = TRUE, scale = 'none', nPcs = 5, seed = 3455)
rpcaDN_seed2 <- pca(t(mat_log_DN_seed2), method = "ppca", 
    center = TRUE, scale = 'none', nPcs = 5, seed = 3455, maxIterations = 1500)
rpcaDMN_seed2 <- pca(t(mat_log_DMN_seed2), method = "ppca", 
    center = TRUE, scale = 'none', nPcs = 5, seed = 3455, maxIterations = 1500)

rpcaNN_leaf2 <- pca(t(mat_log_leaf2), method = "ppca", 
    center = TRUE, scale = 'none', nPcs = 5, seed = 3455)
rpcaDN_leaf2 <- pca(t(mat_log_DN_leaf2), method = "ppca", 
    center = TRUE, scale = 'none', nPcs = 5, seed = 3455, maxIterations = 1500)
rpcaDMN_leaf2 <- pca(t(mat_log_DMN_leaf2), method = "ppca", 
    center = TRUE, scale = 'none', nPcs = 5, seed = 3455, maxIterations = 1500)

# define colors (experiment dependent)
cls_seed1 <- as.character(met1_norm$smp[, 2])
cls_seed2 <- as.character(met2_seed_norm$smp[, 2])
cls_leaf2 <- as.character(met2_leaf_norm$smp[, 2])

pdf('PCA_Mix_sample_seed1_NN_positive.pdf', 7, 7)
slplot(rpcaNN_seed1, pcs = c(1, 2), scoresLoadings = c(TRUE, FALSE), scol = cls_seed1, scex = 1.2)
slplot(rpcaNN_seed1, pcs = c(1, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed1, scex = 1.2)
slplot(rpcaNN_seed1, pcs = c(2, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed1, scex = 1.2)
dev.off()
pdf('PCA_Mix_sample_seed2_NN_positive.pdf', 7, 7)
slplot(rpcaNN_seed2, pcs = c(1, 2), scoresLoadings = c(TRUE, FALSE), scol = cls_seed2, scex = 1.2)
slplot(rpcaNN_seed2, pcs = c(1, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed2, scex = 1.2)
slplot(rpcaNN_seed2, pcs = c(2, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed2, scex = 1.2)
dev.off()
pdf('PCA_Mix_sample_leaf2_NN_positive.pdf', 7, 7)
slplot(rpcaNN_leaf2, pcs = c(1, 2), scoresLoadings = c(TRUE, FALSE), scol = cls_leaf2, scex = 1.2)
slplot(rpcaNN_leaf2, pcs = c(1, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_leaf2, scex = 1.2)
slplot(rpcaNN_leaf2, pcs = c(2, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_leaf2, scex = 1.2)
dev.off()

pdf('PCA_Mix_sample_seed1_DN_positive.pdf', 7, 7)
slplot(rpcaDN_seed1, pcs = c(1, 2), scoresLoadings = c(TRUE, FALSE), scol = cls_seed1, scex = 1.2)
slplot(rpcaDN_seed1, pcs = c(1, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed1, scex = 1.2)
slplot(rpcaDN_seed1, pcs = c(2, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed1, scex = 1.2)
dev.off()
pdf('PCA_Mix_sample_seed2_DN_positive.pdf', 7, 7)
slplot(rpcaDN_seed2, pcs = c(1, 2), scoresLoadings = c(TRUE, FALSE), scol = cls_seed2, scex = 1.2)
slplot(rpcaDN_seed2, pcs = c(1, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed2, scex = 1.2)
slplot(rpcaDN_seed2, pcs = c(2, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed2, scex = 1.2)
dev.off()
pdf('PCA_Mix_sample_leaf2_DN_positive.pdf', 7, 7)
slplot(rpcaDN_leaf2, pcs = c(1, 2), scoresLoadings = c(TRUE, FALSE), scol = cls_leaf2, scex = 1.2)
slplot(rpcaDN_leaf2, pcs = c(1, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_leaf2, scex = 1.2)
slplot(rpcaDN_leaf2, pcs = c(2, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_leaf2, scex = 1.2)
dev.off()

pdf('PCA_Mix_sample_seed1_DMN_positive.pdf', 7, 7)
slplot(rpcaDMN_seed1, pcs = c(1, 2), scoresLoadings = c(TRUE, FALSE), scol = cls_seed1, scex = 1.2)
slplot(rpcaDMN_seed1, pcs = c(1, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed1, scex = 1.2)
slplot(rpcaDMN_seed1, pcs = c(2, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed1, scex = 1.2)
dev.off()
pdf('PCA_Mix_sample_seed2_DMN_positive.pdf', 7, 7)
slplot(rpcaDMN_seed2, pcs = c(1, 2), scoresLoadings = c(TRUE, FALSE), scol = cls_seed2, scex = 1.2)
slplot(rpcaDMN_seed2, pcs = c(1, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed2, scex = 1.2)
slplot(rpcaDMN_seed2, pcs = c(2, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_seed2, scex = 1.2)
dev.off()
pdf('PCA_Mix_sample_leaf2_DMN_positive.pdf', 7, 7)
slplot(rpcaDMN_leaf2, pcs = c(1, 2), scoresLoadings = c(TRUE, FALSE), scol = cls_leaf2, scex = 1.2)
slplot(rpcaDMN_leaf2, pcs = c(1, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_leaf2, scex = 1.2)
slplot(rpcaDMN_leaf2, pcs = c(2, 3), scoresLoadings = c(TRUE, FALSE), scol = cls_leaf2, scex = 1.2)
dev.off()


# save loadings
write.table(loadings(rpcaNN_seed1), file = 'PCA-1000_Loadings_seed1_NN_positive.txt', sep = "\t", quote = FALSE)
write.table(loadings(rpcaNN_seed2), file = 'PCA-1000_Loadings_seed2_NN_positive.txt', sep = "\t", quote = FALSE)
write.table(loadings(rpcaNN_leaf2), file = 'PCA-1000_Loadings_leaf2_NN_positive.txt', sep = "\t", quote = FALSE)

write.table(loadings(rpcaDN_seed1), file = 'PCA-1000_Loadings_seed1_DN_positive.txt', sep = "\t", quote = FALSE)
write.table(loadings(rpcaDN_seed2), file = 'PCA-1000_Loadings_seed2_DN_positive.txt', sep = "\t", quote = FALSE)
write.table(loadings(rpcaDN_leaf2), file = 'PCA-1000_Loadings_leaf2_DN_positive.txt', sep = "\t", quote = FALSE)

write.table(loadings(rpcaDMN_seed1), file = 'PCA-1000_Loadings_seed1_DMN_positive.txt', sep = "\t", quote = FALSE)
write.table(loadings(rpcaDMN_seed2), file = 'PCA-1000_Loadings_seed2_DMN_positive.txt', sep = "\t", quote = FALSE)
write.table(loadings(rpcaDMN_leaf2), file = 'PCA-1000_Loadings_leaf2_DMN_positive.txt', sep = "\t", quote = FALSE)

################################################################################
## match sample numbers to ecotype.number
namesGWAS <- read.xlsx("match_accession_names_to_id.xlsx", sheetIndex = 1)

rename <- function(mat, met_names, weight, split1 = "_Positive", 
    split2 = "Sample", tissue = "seed") {
    
    ## remove columns that contain QC
    if (length(grep(colnames(mat), pattern = "QC")) > 0) 
        mat <- mat[, -grep(colnames(mat), pattern = "QC")]
    
    ind_cut_map <- unlist(lapply(strsplit(colnames(mat), split = split1), "[", 1))
    ind_cut_map <- unlist(lapply(strsplit(ind_cut_map, split = split2), "[", 2))
    
    if (tissue == "seed") {
        ind_to_genotype_map <- ifelse(ind_cut_map %in% 
            weight[, "Internal.ID"], match(ind_cut_map, weight[, "Internal.ID"]), FALSE)
        namesSamples <- unlist(lapply(strsplit(as.character(
            weight[ind_to_genotype_map, "Name"]), split = "[.]"), "[", 2))
    }
    if (tissue == "leaf") {
        ind_to_genotype_map <- ifelse(ind_cut_map %in% 
            weight[, "NO"], match(ind_cut_map, weight[, "NO"]), FALSE)
        namesSamples <- as.character(weight[ind_to_genotype_map, "Accession.name"])
        ##namesSamples <- unlist(lapply(strsplit(as.character(weight[ind_to_genotype_map, "Accession.name"]), split = "[.]"), "[", 2))
    }
    
    ## rename namesSamples that they match to namesGWAS
    namesSamples[which(namesSamples == "Ll-0")] <- "LL-0"
    namesSamples[which(namesSamples == "11PNA4")] <- "11PNA4.101"
    namesSamples[which(namesSamples == "11ME1")] <- "11ME1.32"
    namesSamples[which(namesSamples == "Jea")] <- "JEA"
    
    mat_new <- mat
    ## rename the columns
    colnames(mat_new)[1:dim(mat_new)[2]] <- as.character(
        namesGWAS[match(namesSamples, namesGWAS[, "name"]), "ecotype.name"])
    
    ## remove the columns with NA values
    mat_new <- mat_new[,-which(is.na(colnames(mat_new)))]
    
    ## transpose and remove first two rows
    mat_new <- t(mat_new)
    mat_new <- cbind(Taxa = rownames(mat_new), mat_new)
    
    colnames(mat_new)[2:ncol(mat_new)] <- as.character(met_names[, "Name"])
    mat_new <- rbind(mz = c("mz", met_names[, "mz"]), mat_new)
    mat_new <- rbind(RT = c("RT", met_names[, "RT"]), mat_new)
    
    return(mat_new)
}

mat_log_DN_seed1_rename <- rename(mat_log_DN_seed1, met_names = met1_dw, 
    weight = weights_seeds1)
mat_log_DMN_seed1_rename <- rename(mat_log_DMN_seed1, met_names = met1_dw, 
    weight = weights_seeds1)

mat_log_DN_seed2_rename <- rename(mat_log_DN_seed2, met_names = met2_dw_rem, 
    weight = weights_seeds2)
mat_log_DMN_seed2_rename <- rename(mat_log_DMN_seed2, met_names = met2_dw_rem, 
    weight = weights_seeds2)

mat_log_DN_leaf2_rename <- rename(mat_log_DN_leaf2, met_names = met2_dw_rem, 
    weight = weights_leaf, tissue = "leaf")
mat_log_DMN_leaf2_rename <- rename(mat_log_DMN_leaf2, met_names = met2_dw_rem,
    weight = weights_leaf, tissue = "leaf")

## Export
write.table(mat_log_DN_seed1_rename, file = 'DN_ecotypes_seed1_positive.txt', 
    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mat_log_DMN_seed1_rename, file = 'DMN_ecotypes_seed1_positive.txt', 
    sep = "\t", quote = FALSE, row.names  =  FALSE)

write.table(mat_log_DN_seed2_rename, file = 'DN_ecotypes_seed2_positive.txt', 
    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mat_log_DMN_seed2_rename, file = 'DMN_ecotypes_seed2_positive.txt', 
    sep = "\t", quote = FALSE, row.names  =  FALSE)

write.table(mat_log_DN_leaf2_rename, file = 'DN_ecotypes_leaf2_positive.txt', 
    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mat_log_DMN_leaf2_rename, file = 'DMN_ecotypes_leaf2_positive.txt', 
    sep = "\t", quote = FALSE, row.names  =  FALSE)
