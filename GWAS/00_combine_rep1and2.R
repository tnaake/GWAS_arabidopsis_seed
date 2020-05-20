## read data
setwd("~/01_GWAS/")

polarity <- "positive"

if (polarity == "negative") {
    rep1 <- read.table("PN_rep1_summedint.gda", sep = "\t", dec = ".", header = TRUE)
    rep2 <- read.table("PN_summed_intensity_rep2_leaf.gda", sep = "\t", dec = ".", header = TRUE)
}
if (polarity == "positive") {
    rep1 <- read.table("PP_rep1_summedint.gda", sep = "\t", dec = ".", header = TRUE)
    rep2 <- read.table("PP_leaf_rep2_summed_intensity.gda", sep = "\t", dec = ".", header = TRUE)
}

## remove the first few lines
rep1 <- rep1[-c(1:3), ]
rep2 <- rep2[-c(1:4), ]

rep1 <- rep1[rep1[, "RT"] > 0.5, ]
rep2 <- rep2[rep2[, "RT"] > 0.5, ]
rep1 <- rep1[rep1[, "RT"] < 16, ]
rep2 <- rep2[rep2[, "RT"] < 16, ]

colnames(rep1)[colnames(rep1) == "m.z"] <- "mz"
colnames(rep2)[colnames(rep2) == "m.z"] <- "mz"

rt.dev.range <- c(-0.30, 0.30)
mz.dev.range <- c(-0.01, 0.01)

res <- vector("list", nrow(rep1))
for(i in 1:nrow(rep1)) {
    mz.diff <- rep2$mz - rep1$mz[i]
    rt.diff <- rep2$RT - rep1$RT[i]
    j <- which(rt.diff > rt.dev.range[1] & rt.diff < rt.dev.range[2] &
                    mz.diff > mz.dev.range[1] & mz.diff < mz.dev.range[2])
    res[[i]] <- if(length(j) == 0) NA else j
}


## obtain values only
if (polarity == "negative") {
    rep1_values <- rep1[, which(colnames(rep1) == "Sample1_Negative"):which(colnames(rep1) == "Sample338_Negative")]
    rep2_values <- rep2[, which(colnames(rep2) == "Sample1_Negative"):which(colnames(rep2) == "Sample338_Negative")]
    colnames(rep1_values) <- unlist(lapply(strsplit(colnames(rep1_values), split = "_Negative"), "[", 1))
    colnames(rep2_values) <- unlist(lapply(strsplit(colnames(rep2_values), split = "_Negative"), "[", 1))
}

if (polarity == "positive") {
    rep1_values <- rep1[, which(colnames(rep1) == "Sample1_Positive"):which(colnames(rep1) == "Sample338_Positive")]
    rep2_values <- rep2[, which(colnames(rep2) == "Sample1_Positive"):which(colnames(rep2) == "Sample338_Positive")]
    rep1_values <- rep1_values[, -which(colnames(rep1_values) == "Sample11_Positive")]
    rep1_values <- rep1_values[, -which(colnames(rep1_values) == "Sample27_Positive")]
    colnames(rep1_values) <- unlist(lapply(strsplit(colnames(rep1_values), split = "_Positive"), "[", 1))
    colnames(rep2_values) <- unlist(lapply(strsplit(colnames(rep2_values), split = "_Positive"), "[", 1))
    rep2_values <- rep2_values[, -which(colnames(rep2_values) == "Sample223")]
}

## convert to matrix
rep1_values <- as.matrix(rep1_values)
mode(rep1_values) <- "numeric"
rep2_values <- as.matrix(rep2_values)
mode(rep2_values) <- "numeric"

## check that we have the same column names
colnames(rep1_values)[!colnames(rep1_values) %in% colnames(rep2_values)]
colnames(rep2_values)[!colnames(rep2_values) %in% colnames(rep1_values)]

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
            covs <- cov(t(log2(rep1_values[rep(i, length(inds)), ] + 1)), t(log2(rep2_values[inds, ] + 1)))
            covs <- covs[1, ]
        } else {
            covs <- cov(log2(rep1_values[i, ] + 1), log2(rep2_values[inds,] + 1))
        }
        cov_max <- covs[which.max(covs)]
        cov_res[i] <- cov_max
        res_cor[[i]] <- res[[i]][which.max(covs)]
    } else {
        cov_res[i] <- NaN
        res_cor[[i]] <- NaN
    }
    
}

res_cor <- unlist(res_cor)
## set features whose features is equal or smaller than 0 to NaN
##res_cor[cov_res <= 0] <- NaN
##cov_res[cov_res <= 0] <- NaN

rep1_cor <- rep1[!is.na(res_cor), ]
rep2_cor <- rep2[res_cor[!is.na(res_cor)], ]

## add the relation of the mapping
rep1_cor <- cbind(rep1_cor, 
    mapping_rep1_rep2 = paste(rep1_cor[, "Name"], rep2[res_cor[!is.na(res_cor)], "Name"], sep = "/"))
rep2_cor <- cbind(rep2_cor, 
    mapping_rep1_rep2 = paste(rep1[!is.na(res_cor), "Name"], rep2_cor[, "Name"], sep = "/"))

## add the covariance with the mapped feature
rep1_cor <- cbind(rep1_cor, cov = cov_res[!is.na(cov_res)])
rep2_cor <- cbind(rep2_cor, cov = cov_res[!is.na(cov_res)])

## add RT deviance 
rep1_cor <- cbind(rep1_cor, rt_dev = rep1_cor[, "RT"] - rep2_cor[, "RT"])
rep2_cor <- cbind(rep2_cor, rt_dev = rep1_cor[, "RT"] - rep2_cor[, "RT"])

## add m/z deviance
rep1_cor <- cbind(rep1_cor, mz_dev = rep1_cor[, "mz"] - rep2_cor[, "mz"])
rep2_cor <- cbind(rep2_cor, mz_dev = rep1_cor[, "mz"] - rep2_cor[, "mz"])


## check if major glucosides and flavonoids are present (use rep2)
check_met <- function(met = rep2_cor, mz_1, mz_2, rt_1, rt_2) {
    rt <- met[, "RT"]
    mz <- met[, "mz"]
    return(met[mz > mz_1 & mz < mz_2 & rt > rt_1 & rt < rt_2, ])
}

add_features <- function(mz_l = 593.13, mz_h = 593.17, rt1_l = 6.40, 
        rt1_h = 6.45, rt2_l = 6.66, rt2_h = 6.70) {
    
    add_rep2 <- rep2[rep2[, "mz"] > mz_l & rep2[, "mz"] < mz_h & 
        rep2[, "RT"] > rt2_l & rep2[, "RT"] < rt2_h, ]
    add_rep1 <- rep1[rep1[, "mz"] > mz_l & rep1[, "mz"] < mz_h & 
        rep1[, "RT"] > rt1_l & rep1[, "RT"] < rt1_h, ]
    
    relation_add <- paste(add_rep1[, "Name"], add_rep2[, "Name"], sep = "/")
    rt_add <- add_rep1[, "RT"] - add_rep2[, "RT"]
    mz_add <- add_rep1[, "mz"] - add_rep2[, "mz"]
    
    d1 <- dim(add_rep1)[1]
    if (d1 > 1) {
        covs <- cov(t(log2(rep2_values[rep(rownames(add_rep2), d1), ] + 1)), 
            t(log2(rep1_values[rownames(add_rep1), ] + 1)))
        covs <- covs[1, ]
    } else {
        covs <- cov(log2(rep2_values[rownames(add_rep2), ] + 1), 
            log2(rep1_values[rownames(add_rep1), ] + 1))
    }
    
    add_rep1_fi <- cbind(add_rep1, mapping_rep1_rep2 = relation_add, cov = covs, 
        rt_dev = rt_add, mz_dev = mz_add)
    add_rep2_fi <- cbind(add_rep2[rep(1, d1), ], mapping_rep1_rep2 = relation_add, 
        cov = covs, rt_dev = rt_add, mz_dev = mz_add)
    
    rownames(add_rep1_fi) <- paste(rownames(add_rep1), "manual", sep = ".")
    rownames(add_rep2_fi) <- paste(rep(rownames(add_rep2), d1), 1:d1, "manual", sep = ".")
    
    return(list(add_rep1_fi, add_rep2_fi))
}

## found
if (polarity == "negative") {
    check_met(mz_1 = 436.00, mz_2 = 436.10, rt_1 = 1.53, rt_2 = 1.57) ## 4-methylsulfinylbutyl glucosinolate (neg)
    check_met(mz_1 = 450.05, mz_2 = 450.10, rt_1 = 2.65, rt_2 = 2.69) ## 5-methylsulfinylpentyl glucosinolate (neg)
    check_met(mz_1 = 464.02, mz_2 = 464.14, rt_1 = 3.36, rt_2 = 3.40) ## 6-methylsulfinylhexyl glucosinolate (neg)
    check_met(mz_1 = 478.06, mz_2 = 478.12, rt_1 = 4.12, rt_2 = 4.16) ## 7-methylsulfinylheptyl glucosinolate (neg)
    check_met(mz_1 = 492.07, mz_2 = 492.13, rt_1 = 4.96, rt_2 = 5.00) ## 8-Methylsulfinyloctyl glucosinolate (neg)
    check_met(mz_1 = 376.02, mz_2 = 376.06, rt_1 = 0.81, rt_2 = 0.87) ## 3-hydroxypropyl-glucosinolate (neg)
    check_met(mz_1 = 390.04, mz_2 = 390.08, rt_1 = 1.18, rt_2 = 1.22) ## 4-hydroxybutylglucosinolate (neg)
    check_met(mz_1 = 494.05, mz_2 = 494.10, rt_1 = 7.26, rt_2 = 7.33) ## 4-benzoyloxybutylglucosinolate (neg)
    check_met(mz_1 = 577.10, mz_2 = 577.20, rt_1 = 7.25, rt_2 = 7.29) ## Kaempferol 3-O-rhamnoside 7-O-rhamnoside
    check_met(mz_1 = 739.21, mz_2 = 739.28, rt_1 = 5.83, rt_2 = 5.87) ## Kaempferol 3-O-[2''-O-(rhamnosyl) glucoside] 7-O-rhamnoside
    check_met(mz_1 = 755.18, mz_2 = 755.28, rt_1 = 5.53, rt_2 = 5.57) ## Quercetin 3-O-[2''-O-(rhamnosyl) glucoside] 7-O-rhamnoside
    check_met(mz_1 = 406.01, mz_2 = 406.05, rt_1 = 3.39, rt_2 = 3.46) ## 3-methylthiopropyl-glucosinolate (neg)
    check_met(mz_1 = 420.02, mz_2 = 420.07, rt_1 = 4.24, rt_2 = 4.28) ## 4-methylthiobutyl glucosinolate (neg)
    check_met(mz_1 = 434.04, mz_2 = 434.09, rt_1 = 5.30, rt_2 = 5.34) ## 5-methylthiopentylglucosinolate (neg)
    check_met(mz_1 = 448.03, mz_2 = 448.10, rt_1 = 6.10, rt_2 = 6.35) ## 6-methylthiohexylglucosinolate (neg)
    check_met(mz_1 = 462.06, mz_2 = 462.12, rt_1 = 7.60, rt_2 = 7.72) ## 7-methylthioheptyl glucosinolate (neg)
    check_met(mz_1 = 476.08, mz_2 = 476.12, rt_1 = 8.82, rt_2 = 8.88) ## 8-Methylthiooctyl glucosinolate (neg)
    check_met(mz_1 = 477.03, mz_2 = 477.10, rt_1 = 5.59, rt_2 = 5.63) ## 4-methoxy-3-indolylmethyl-glucosinolate (neg)
    check_met(mz_1 = 480.03, mz_2 = 480.09, rt_1 = 6.43, rt_2 = 6.47) ## 3-benzoyloxypropyl-glucosinolate (neg)
    check_met(mz_1 = 477.03, mz_2 = 477.10, rt_1 = 6.49, rt_2 = 6.55) ## 1-methoxy-3-indolylmethyl-glucosinolate (neg)
    check_met(mz_1 = 582.06, mz_2 = 582.13, rt_1 = 6.62, rt_2 = 6.66) ## 3-sinapoyloxypropylglucosinolate (neg)

    ## not found
    check_met(mz_1 = 593.13, mz_2 = 593.17, rt_1 = 6.66, rt_2 = 6.70) ## Kaempferol 3-O-glucoside 7-O-rhamnoside
    check_met(mz_1 = 593.13, mz_2 = 593.17, rt_1 = 6.75, rt_2 = 6.79) ## Quercetin 3-O-rhamnoside 7-O-rhamnoside
    check_met(mz_1 = 755.18, mz_2 = 755.24, rt_1 = 6.44, rt_2 = 6.48) ## Kaempferol 3-O-glucosyl-glucoside 7-O-rhamnoside
    check_met(mz_1 = 609.12, mz_2 = 609.20, rt_1 = 6.15, rt_2 = 6.30) ## Quercetin 3-O-glucoside 7-O-rhamnoside

    ## Kaempferol 3-O-glucoside 7-O-rhamnoside
    rep1_cor <- rbind(rep1_cor, add_features(mz_l = 593.13, mz_h = 593.17, rt1_l = 6.40, rt1_h = 6.45, rt2_l = 6.66, rt2_h = 6.70)[[1]])
    rep2_cor <- rbind(rep2_cor, add_features(mz_l = 593.13, mz_h = 593.17, rt1_l = 6.40, rt1_h = 6.45, rt2_l = 6.66, rt2_h = 6.70)[[2]])
    ## Quercetin 3-O-rhamnoside 7-O-rhamnoside
    rep1_cor <- rbind(rep1_cor, add_features(mz_l = 593.13, mz_h = 593.17, rt1_l = 6.60, rt1_h = 6.65, rt2_l = 6.75, rt2_h = 6.79)[[1]])
    rep2_cor <- rbind(rep2_cor, add_features(mz_l = 593.13, mz_h = 593.17, rt1_l = 6.60, rt1_h = 6.65, rt2_l = 6.75, rt2_h = 6.79)[[2]])
    ## Kaempferol 3-O-glucosyl-glucoside 7-O-rhamnoside
    rep1_cor <- rbind(rep1_cor, add_features(mz_l = 755.18, mz_h = 755.23, rt1_l = 6.10, rt1_h = 6.40, rt2_l = 6.44, rt2_h = 6.48)[[1]])
    rep2_cor <- rbind(rep2_cor, add_features(mz_l = 755.18, mz_h = 755.23, rt1_l = 6.10, rt1_h = 6.40, rt2_l = 6.44, rt2_h = 6.48)[[2]])
    ## Quercetin 3-O-glucoside 7-O-rhamnoside
    rep1_cor <- rbind(rep1_cor, add_features(mz_l = 609.12, mz_h = 609.20, rt1_l = 6.16, rt1_h = 6.20, rt2_l = 6.20, rt2_h = 6.26)[[1]])
    rep2_cor <- rbind(rep2_cor, add_features(mz_l = 609.12, mz_h = 609.20, rt1_l = 6.16, rt1_h = 6.20, rt2_l = 6.20, rt2_h = 6.26)[[2]])
}

if (polarity == "positive") {
    check_met(mz_1 = 476.00, mz_2 = 476.05, rt_1 = 1.53, rt_2 = 1.60) ## 4-methylsulfinylbutyl glucosinolate (pos)
    check_met(mz_1 = 490.00, mz_2 = 490.08, rt_1 = 2.65, rt_2 = 2.69) ## 5-methylsulfinylpentyl glucosinolate (pos)
    check_met(mz_1 = 504.00, mz_2 = 504.10, rt_1 = 3.36, rt_2 = 3.40) ## 6-methylsulfinylhexyl glucosinolate (pos)
    check_met(mz_1 = 518.00, mz_2 = 518.11, rt_1 = 4.35, rt_2 = 4.46) ## 7-methylsulfinylheptyl glucosinolate (pos)
    check_met(mz_1 = 532.02, mz_2 = 532.13, rt_1 = 5.60, rt_2 = 5.70) ## 8-Methylsulfinyloctyl glucosinolate (pos)
    check_met(mz_1 = 415.98, mz_2 = 416.02, rt_1 = 0.81, rt_2 = 0.87) ## 3-hydroxypropyl-glucosinolate (pos)
    check_met(mz_1 = 430.00, mz_2 = 430.06, rt_1 = 1.18, rt_2 = 1.22) ## 4-hydroxybutylglucosinolate (pos)
    check_met(mz_1 = 424.08, mz_2 = 424.13, rt_1 = 6.43, rt_2 = 6.6) ## 3-benzoyloxypropyl-glucosinolate (pos)
    check_met(mz_1 = 534.02, mz_2 = 534.10, rt_1 = 7.26, rt_2 = 7.33) ## 4-benzoyloxybutylglucosinolate (pos)
    check_met(mz_1 = 579.16, mz_2 = 579.20, rt_1 = 7.25, rt_2 = 7.29) ## Kaempferol 3-O-rhamnoside 7-O-rhamnoside
    check_met(mz_1 = 741.20, mz_2 = 741.25, rt_1 = 5.83, rt_2 = 5.87) ## Kaempferol 3-O-[2''-O-(rhamnosyl) glucoside] 7-O-rhamnoside
    check_met(mz_1 = 757.19, mz_2 = 757.25, rt_1 = 5.40, rt_2 = 5.50) ## Quercetin 3-O-[2''-O-(rhamnosyl) glucoside] 7-O-rhamnoside
    check_met(mz_1 = 445.98, mz_2 = 446.05, rt_1 = 3.39, rt_2 = 3.46) ## 3-methylthiopropyl-glucosinolate (pos)
    check_met(mz_1 = 460.00, mz_2 = 460.04, rt_1 = 4.24, rt_2 = 4.28) ## 4-methylthiobutyl glucosinolate (pos)
    check_met(mz_1 = 378.08, mz_2 = 378.15, rt_1 = 5.25, rt_2 = 5.34) ## 5-methylthiopentylglucosinolate (pos)
    check_met(mz_1 = 488.01, mz_2 = 488.10, rt_1 = 6.20, rt_2 = 6.30) ## 6-methylthiohexylglucosinolate (pos)
    check_met(mz_1 = 406.10, mz_2 = 406.17, rt_1 = 7.60, rt_2 = 7.72) ## 7-methylthioheptyl glucosinolate (pos)
    check_met(mz_1 = 420.12, mz_2 = 420.19, rt_1 = 8.82, rt_2 = 8.88) ## 8-Methylthiooctyl glucosinolate (pos)
    check_met(mz_1 = 479.05, mz_2 = 479.09, rt_1 = 5.59, rt_2 = 5.63) ## 4-methoxy-3-indolylmethyl-glucosinolate (pos)
    check_met(mz_1 = 595.13, mz_2 = 595.19, rt_1 = 6.46, rt_2 = 6.70) ## Kaempferol 3-O-glucoside 7-O-rhamnoside
    check_met(mz_1 = 595.13, mz_2 = 595.19, rt_1 = 6.75, rt_2 = 6.79) ## Quercetin 3-O-rhamnoside 7-O-rhamnoside
    check_met(mz_1 = 757.18, mz_2 = 757.24, rt_1 = 6.44, rt_2 = 6.48) ## Kaempferol 3-O-glucosyl-glucoside 7-O-rhamnoside
    
    ## not found
    check_met(mz_1 = 479.05, mz_2 = 479.09, rt_1 = 6.49, rt_2 = 6.55) ## 1-methoxy-3-indolylmethyl-glucosinolate (pos)
    check_met(mz_1 = 611.12, mz_2 = 611.20, rt_1 = 6.15, rt_2 = 6.30) ## Quercetin 3-O-glucoside 7-O-rhamnoside
    
    ## 5-methylthiopentylglucosinolate (pos)
    rep1_cor <- rbind(rep1_cor, add_features(mz_l = 479.05, mz_h = 479.09, rt1_l = 6.80, rt1_h = 6.90, rt2_l = 6.43, rt2_h = 6.55)[[1]])
    rep2_cor <- rbind(rep2_cor, add_features(mz_l = 479.05, mz_h = 479.09, rt1_l = 6.80, rt1_h = 6.90, rt2_l = 6.43, rt2_h = 6.55)[[2]])
    
    ## Quercetin 3-O-glucoside 7-O-rhamnoside
    rep1_cor <- rbind(rep1_cor, add_features(mz_l = 611.12, mz_h = 611.20, rt1_l = 6.00, rt1_h = 6.20, rt2_l = 6.20, rt2_h = 6.30)[[1]])
    rep2_cor <- rbind(rep2_cor, add_features(mz_l = 611.12, mz_h = 611.20, rt1_l = 6.00, rt1_h = 6.20, rt2_l = 6.20, rt2_h = 6.30)[[2]])
    
}


## save
write.table(rep1_cor, file = paste0("./rep1_", polarity, "_match.csv"), sep = "\t", dec = ".", quote = FALSE)
write.table(rep2_cor, file = paste0("./rep2_", polarity, "_match.csv"), sep = "\t", dec = ".", quote = FALSE)

