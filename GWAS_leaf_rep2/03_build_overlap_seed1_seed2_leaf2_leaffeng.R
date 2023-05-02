## change path
options(stringsAsFactors = FALSE)
setwd("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/")

## load the file from Feng
leaf_feng_normalized <- read.delim(
    "files_feng/combine_results/normalized/normalized_gwas_gene_info.txt", 
    header = TRUE, sep = "\t", dec = ".", quote = "") ## dim 17835 210 ## 36107 47
leaf_feng_normalized_loci_info <- read.delim(
    "files_feng/combine_results/normalized/normalized_gwas_loci_info.txt",
    header = TRUE,sep = "\t", dec = ".", quote = "")
leaf_feng_batch <- read.delim(
    "files_feng/combine_results/normalized_batch_corrected/normalized_batch_corrected_gwas_gene_info.txt", 
    header = TRUE, sep = "\t", dec = ".", quote = "") ## dim 17553 152 ## 36107 47
leaf_feng_batch_loci_info <- read.delim(
    "files_feng/combine_results/normalized_batch_corrected/normalized_batch_corrected_gwas_locus_info.txt", 
    header = TRUE, sep = "\t", dec = ".", quote = "") ## dim 17553 152 ## 36107 47


## load the gwas_complete_met_all_trueLociLOD_neg.txt file that contains the
## relations between seed1, seed2, leaf2
trueLociLOD <- read.delim("../gwas_complete_met_all_trueLociLOD_neg.txt",
    header = TRUE, sep = "\t", dec = ".", quote = "")

## load the relation to rep1, rep2, rep_leaf: there is a link between rep2 and 
## rep_feng (column mapping_rep2_rep_feng)
rel_rep1 <- read.csv("rep1_negative_match.csv", sep = "\t")
rel_rep2 <- read.csv("rep2_negative_match.csv", sep = "\t")
rel_rep_feng <- read.csv("rep_feng_negative_match.csv", sep = "\t")

## obtain the mapping relation
mapping_rep2 <- lapply(
    strsplit(rel_rep_feng[["mapping_rep2_rep_feng"]], split = "/"), "[", 1) |>
    unlist()
names(mapping_rep2) <- lapply(
    strsplit(rel_rep_feng[["mapping_rep2_rep_feng"]], split = "/"), "[", 2) |>
    unlist()
mapping_rep2 <- mapping_rep2[
    !duplicated(paste(names(mapping_rep2), as.character(mapping_rep2)))]


add_to_trueLociLOD <- function(gene_info = leaf_feng, loci_info = leaf_feng_loci_info, trueLociLOD = trueLociLOD) {
    
    ## 1) add columns to trueLociLOD
    trueLociLOD[["locusID_leaf_feng"]] <- NA
    trueLociLOD[["bestSNP_lod_leaf_feng"]] <- NA
    trueLociLOD[["locus_tag_leaf_feng"]] <- NA
    trueLociLOD[["met_leaf_feng"]] <- NA
    trueLociLOD[["mz_leaf_feng"]] <- NA
    trueLociLOD[["rt_leaf_feng"]] <- NA
    
    ## 2) create object to store results
    trueLociLOD_res <- trueLociLOD
    
    ## 3) find columns with mass features, remove the "MLM." in front of 
    ## "Cluster"
    cols_features <- grep(colnames(gene_info), pattern = "Peak[.]ID_")
    gene_info[, cols_features] <- lapply(gene_info[, cols_features], 
        gsub, pattern = "MLM[.]", replacement = "")
    
    ## 4) remove the column Function, product, and note
    gene_info <- gene_info[, !colnames(gene_info) %in% 
        c("Function", "product", "note")]
    cols_features <- grep(colnames(gene_info), pattern = "Peak[.]ID_")
    
    ## 5) create three lists from trueLociLOD that contain the loci regions:
    ## each entry contains the chromosome, the start and end gene
    get_loci <- function(loci) {
        lapply(loci, function(i) {
            if (length(i) != 0) {
                region <- strsplit(i, split = "AT") |> 
                    unlist()
                region <- region[region != ""]
                strsplit(region, split = "G") |>
                    unlist() |>
                    as.numeric()
            } else NaN
        })   
    }
    loci_seed1 <- strsplit(trueLociLOD[["locus_tag_seed1"]], split = "[.][.]")
    loci_seed1 <- get_loci(loci_seed1)
    loci_seed1_chr <- unlist(lapply(loci_seed1, "[", 1))
    loci_seed1_min <- unlist(lapply(loci_seed1, "[", 2))
    loci_seed1_max <- unlist(lapply(loci_seed1, "[", 4))
    loci_seed2 <- strsplit(trueLociLOD[["locus_tag_seed2"]], split = "[.][.]")
    loci_seed2 <- get_loci(loci_seed2)
    loci_seed2_chr <- unlist(lapply(loci_seed2, "[", 1))
    loci_seed2_min <- unlist(lapply(loci_seed2, "[", 2))
    loci_seed2_max <- unlist(lapply(loci_seed2, "[", 4))
    loci_leaf2 <- strsplit(trueLociLOD[["locus_tag_leaf2"]], split = "[.][.]")
    loci_leaf2 <- get_loci(loci_leaf2)
    loci_leaf2_chr <- unlist(lapply(loci_leaf2, "[", 1))
    loci_leaf2_min <- unlist(lapply(loci_leaf2, "[", 2))
    loci_leaf2_max <- unlist(lapply(loci_leaf2, "[", 4))
    
    ## 6) iterate through unique loci in gene_info
    ## 6.1) create unique loci
    gene_info$unique_loci <- paste(gene_info$locusID, gene_info$no_SNP, 
        gene_info$best_SNP, gene_info$chrom, gene_info$best_SNP_lod, sep = "_")
    
    ## 6.2) find unique loci
    unique_loci_u <- gene_info$unique_loci |>
        unique()
    
    ## 6.3) iterate through unique loci and check overlap with trueLociLOD
    for (i in unique_loci_u) {
        
        gene_info_loci <- gene_info[gene_info$unique_loci == i, ]
        
        ## find mass features
        features_i <- gene_info_loci[, cols_features] |> 
            unlist() |>
            unique()
        features_i <- features_i[features_i != "."]
        
        ## get the positions of the locus tag
        gene_pos_i <- lapply(get_loci(gene_info_loci[["locus_tag"]]), "[", 2) |>
            unlist()
        
        ## check if the loci are overlapping: this will be the case if 
        ## either the min of loci_* lies within gene_pos_i OR the max of
        ## loci_* lies within gene_pos_i OR if gene_pos_i lies within loci_*
        ## do this separately for seed1, seed2, leaf2
        ind_loci_seed1 <- (
            (loci_seed1_min <= max(gene_pos_i) & loci_seed1_min >= min(gene_pos_i)) |
                (loci_seed1_max >= min(gene_pos_i) & loci_seed1_max <= max(gene_pos_i)) |
                (loci_seed1_min <= min(gene_pos_i) & loci_seed1_max >= max(gene_pos_i))
        ) & loci_seed1_chr == gene_info_loci[["chrom"]][1]
        ind_loci_seed2 <- (
            (loci_seed2_min <= max(gene_pos_i) & loci_seed2_min >= min(gene_pos_i)) |
                (loci_seed2_max >= min(gene_pos_i) & loci_seed2_max <= max(gene_pos_i)) |
                (loci_seed2_min <= min(gene_pos_i) & loci_seed2_max >= max(gene_pos_i))
        ) & loci_seed2_chr == gene_info_loci[["chrom"]][1]
        ind_loci_leaf2 <- (
            (loci_leaf2_min <= max(gene_pos_i) & loci_leaf2_min >= min(gene_pos_i)) |
                (loci_leaf2_max >= min(gene_pos_i) & loci_leaf2_max <= max(gene_pos_i)) |
                (loci_leaf2_min <= min(gene_pos_i) & loci_leaf2_max >= max(gene_pos_i))
        ) & loci_leaf2_chr == gene_info_loci[["chrom"]][1]
        
        ## iterate through the features_i
        for (features_i_j in features_i) {
            
            ## get name of feature in rep2
            ind_feature <- trueLociLOD[["met_rep2"]] == mapping_rep2[features_i_j]
            
            ## find the rows that both have the sampe feature rep2 and 
            ## overlapping loci
            ind <- ind_feature & (ind_loci_seed1 | ind_loci_seed2 | ind_loci_leaf2)
            ind <- which(ind)
            
            ## if no overlap, write to end
            if (length(ind) == 0) {
                ind <- nrow(trueLociLOD_res) + 1
                trueLociLOD_res[ind, ] <- rep(NA, 21)
            }
            
            ## write the information to trueLociLOD_res
            trueLociLOD_res[["locusID_leaf_feng"]][ind] <- 
                unique(gene_info_loci[["unique_loci"]])
            trueLociLOD_res[["bestSNP_lod_leaf_feng"]][ind] <-
                max(loci_info[
                    loci_info[["locusID"]] %in% unique(gene_info_loci[["locusID"]]) &
                    loci_info[["Peak.ID"]] %in% paste0("MLM.", features_i_j), "lod"], 
                na.rm = TRUE)
                ##unique(gene_info_loci[["best_SNP_lod"]])
            trueLociLOD_res[["locus_tag_leaf_feng"]][ind] <-
                paste(gene_info_loci[["locus_tag"]][1], 
                    gene_info_loci[["locus_tag"]][nrow(gene_info_loci)], 
                    sep = "..")
            trueLociLOD_res[["met_leaf_feng"]][ind] <- features_i_j
            trueLociLOD_res[["mz_leaf_feng"]][ind] <- 
                unique(rel_rep_feng[rel_rep_feng[["Name"]] == features_i_j, "mz"])
            trueLociLOD_res[["rt_leaf_feng"]][ind] <-
                unique(rel_rep_feng[rel_rep_feng[["Name"]] == features_i_j, "RT"])
        }
    }
    return(trueLociLOD_res)
}

## apply the function on leaf_feng_normalized and leaf_feng_batch
trueLociLOD_normalized <- add_to_trueLociLOD(gene_info = leaf_feng_normalized, 
    loci_info = leaf_feng_normalized_loci_info, trueLociLOD = trueLociLOD)
trueLociLOD_batch <- add_to_trueLociLOD(gene_info = leaf_feng_batch, 
    loci_info = leaf_feng_batch_loci_info, trueLociLOD = trueLociLOD)

write.table(trueLociLOD_normalized, 
    file = "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_rep12_normalized_neg.txt",
    sep = "\t", dec = ".", quote = FALSE)
write.table(trueLociLOD_batch,
    file = "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_rep12_normalized_batch_corrected.txt", 
    sep = "\t", dec = ".", quote = FALSE)
