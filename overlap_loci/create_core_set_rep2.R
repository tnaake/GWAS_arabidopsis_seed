## negative mode

## load trueLociLOD files 
trueLociLOD <- read.table( 
    file = "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_rep12_normalized_neg.txt",
    header = TRUE, sep = "\t", dec = ".", quote = "")
##trueLociLOD <- read.delim("../gwas_complete_met_all_trueLociLOD_neg.txt",
##    header = TRUE, sep = "\t", dec = ".", quote = "")


## load the relation to rep1, rep2, rep_leaf: there is a link between rep2 and 
## rep_feng (column mapping_rep2_rep_feng)
rel_rep1 <- read.csv("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/rep1_negative_match.csv", sep = "\t")
rel_rep_feng <- read.csv("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/rep_feng_negative_match.csv", sep = "\t")

## determine the covariance that results in 9008 features after the truncation step
features <- rel_rep1[rel_rep1$cov > 1.3874,  "mapping_rep1_rep2"]
features_rep1 <- lapply(strsplit(features, split = "/"), "[", 1) |>
    unlist()
features_rep2 <- lapply(strsplit(features, split = "/"), "[", 2) |>
    unlist()

## what are those features in features_feng?
features_rep_feng_rep2 <- lapply(strsplit(rel_rep_feng[["mapping_rep2_rep_feng"]], split = "/"), "[", 1) |>
    unlist()
features_rep_feng_rep_feng <- lapply(strsplit(rel_rep_feng[["mapping_rep2_rep_feng"]], split = "/"), "[", 2) |>
    unlist()
features_rep_feng <- features_rep_feng_rep_feng[features_rep_feng_rep2 %in% features_rep2]


## truncate
trueLociLOD_cut <- trueLociLOD |>
    dplyr::filter(met_rep1 %in% features_rep1 | is.na(met_rep1)) |>
    dplyr::filter(met_rep2 %in% features_rep2 | is.na(met_rep2)) |>
    dplyr::filter(met_leaf_feng %in% features_rep_feng | is.na(met_leaf_feng))


## write to file
write.table(trueLociLOD_cut, 
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
features <- rel_rep1[rel_rep1$cov > 2.8474, "mapping_rep1_rep2"]
features_rep1 <- lapply(strsplit(features, split = "/"), "[", 1) |>
    unlist()
features_rep2 <- lapply(strsplit(features, split = "/"), "[", 2) |>
    unlist()

## truncate
trueLociLOD_cut <- trueLociLOD |>
    dplyr::filter(met_rep1 %in% features_rep1) |>
    dplyr::filter(met_rep2 %in% features_rep2)

## write to file
write.table(trueLociLOD_cut, 
    "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_pos_cut.txt", 
    sep = "\t", quote = FALSE)
