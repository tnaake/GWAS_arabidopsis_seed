## negative mode

## load trueLociLOD files 
trueLociLOD_normalized <- read.table( 
    file = "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_rep12_normalized_neg.txt",
    header = TRUE, sep = "\t", dec = ".", quote = "")
##trueLociLOD <- read.delim("../gwas_complete_met_all_trueLociLOD_neg.txt",
##    header = TRUE, sep = "\t", dec = ".", quote = "")


## load the relation to rep1, rep2, rep_leaf: there is a link between rep2 and 
## rep_feng (column mapping_rep2_rep_feng)
rel_rep1 <- read.csv("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/rep1_negative_match.csv", sep = "\t")

## determine the covariance that results in 9008 features
features <- rel_rep1[rel_rep1$cov > 1.5357, "mapping_rep1_rep2"]
features_rep1 <- lapply(strsplit(features, split = "/"), "[", 1) |>
    unlist()
features_rep2 <- lapply(strsplit(features, split = "/"), "[", 2) |>
    unlist()

## truncate
trueLociLOD_normalized_cut <- trueLociLOD_normalized[
    (trueLociLOD_normalized$met_rep1 %in% features_rep1 &
         trueLociLOD_normalized$met_rep2 %in% features_rep2) |
        (is.na(trueLociLOD_normalized$met_rep1) & 
            is.na(trueLociLOD_normalized$met_rep2)), ]

## write to file
write.table(trueLociLOD_normalized_cut, 
    "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_neg_cut.txt", 
    sep = "\t", quote = FALSE)

## positive mode

## load trueLociLOD files 
trueLociLOD <- read.delim(
    "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_pos.txt",
    header = TRUE, sep = "\t", dec = ".", quote = "")

## load the relation to rep1, rep2, rep_leaf: there is a link between rep2 and 
## rep_feng (column mapping_rep2_rep_feng)
rel_rep1 <- read.csv("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/rep1_positive_match.csv", sep = "\t")

## determine the covariance that results in 12133 features
features <- rel_rep1[rel_rep1$cov > 3.168, "mapping_rep1_rep2"]
features_rep1 <- lapply(strsplit(features, split = "/"), "[", 1) |>
    unlist()
features_rep2 <- lapply(strsplit(features, split = "/"), "[", 2) |>
    unlist()

## truncate
trueLociLOD_cut <- trueLociLOD[
    trueLociLOD$met_rep1 %in% features_rep1 & 
        trueLociLOD$met_rep2 %in% features_rep2, ]

## write to file
write.table(trueLociLOD_cut, 
    "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_pos_cut.txt", 
    sep = "\t", quote = FALSE)
