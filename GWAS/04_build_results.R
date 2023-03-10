## change path
options(stringsAsFactors = FALSE)
setwd("~/GitHub/GWAS_arabidopsis_seed/results_build_loci")

polarity <- "pos"
## change 
if (polarity == "neg") {
    seed1 <- read.delim("gwas_gene_info_seed_1_neg_mod.txt", header = TRUE, 
        sep = "\t", dec = ".") ## dim 32434  4533
    seed2 <- read.delim("gwas_gene_info_seed_2_neg_mod.txt", header = TRUE, 
        sep = "\t", dec = ".") ## dim 32079 2274
    leaf2 <- read.delim("gwas_gene_info_leaf_neg_mod.txt", header = TRUE, 
        sep = "\t", dec = ".") ## dim 32199 1905
}

if (polarity == "pos") {
    seed1 <- read.delim("gwas_gene_info_seed_1_pos_mod.txt", header = TRUE, 
        sep = "\t", dec = ".") ## dim 32486 7000
    seed2 <- read.delim("gwas_gene_info_seed_2_pos_mod.txt", header = TRUE, 
        sep = "\t", dec = ".") ## dim 32331 4245
    leaf2 <- read.delim("gwas_gene_info_leaf_pos_mod.txt", header = TRUE, 
        sep = "\t", dec = ".") ## dim 32344 4090
}

## load file with relation about metabolites (contains also information on m/z and rt of each metabolite)
setwd("~/GitHub/GWAS_arabidopsis_seed/")

if (polarity == "neg") {
    rel_rep1 <- read.csv("rep1_negative_match.csv", sep = "\t")
    rel_rep2 <- read.csv("rep2_negative_match.csv", sep = "\t")
} 
if (polarity == "pos") {
    rel_rep1 <- read.csv("rep1_positive_match.csv", sep = "\t")
    rel_rep2 <- read.csv("rep2_positive_match.csv", sep = "\t")
}
all(rel_rep1[, "mapping_rep1_rep2"] == rel_rep2[, "mapping_rep1_rep2"])

## retrieve the information
relat <- rel_rep1[, "mapping_rep1_rep2"]
relat <- strsplit(relat, split = "/")
relat <- cbind(unlist(lapply(relat, "[", 1)), 
    unlist(lapply(relat, "[", 2)), unlist(lapply(relat, "[", 2)))
colnames(relat) <- c("seed1", "seed2", "leaf")

## add MLM. to beginning
relat[, 1] <- paste("MLM.", relat[, 1], sep = "")
relat[, 2] <- paste("MLM.", relat[, 2], sep = "")
relat[, 3] <- paste("MLM.", relat[, 3], sep = "")

## get all mapped metabolites
seed1_peak <- seed1[, grep(colnames(seed1), pattern = "Peak")]
seed1_peak <- unlist(seed1_peak)
seed1_ch <- as.character(seed1_peak)

seed2_peak <- seed2[, grep(colnames(seed2), pattern = "Peak")]
seed2_peak <- unlist(seed2_peak)
seed2_ch <- as.character(seed2_peak)

leaf2_peak <- leaf2[, grep(colnames(leaf2), pattern = "Peak")]
leaf2_peak <- unlist(leaf2_peak)
leaf2_ch <- as.character(leaf2_peak)

## make unique
metab_seed1_uniq <- seed1_ch[grep(seed1_ch, pattern = "MLM")]
metab_seed1_uniq <- as.character(metab_seed1_uniq)
metab_seed1_uniq <- unique(metab_seed1_uniq)

metab_seed2_uniq <- seed2_ch[grep(seed2_ch, pattern = "MLM")]
metab_seed2_uniq <- as.character(metab_seed2_uniq)
metab_seed2_uniq <- unique(metab_seed2_uniq)

metab_leaf2_uniq <- leaf2_ch[grep(leaf2_ch, pattern = "MLM")]
metab_leaf2_uniq <- as.character(metab_leaf2_uniq)
metab_leaf2_uniq <- unique(metab_leaf2_uniq)

met_all <- vector("list", nrow(relat))


## assume metab_seed1_uniq[1], metab_seed2_uniq[1], metab_leaf2_uniq[1] belong together
res <- matrix(0, ncol = 15, nrow = 0)
colnames(res) <- c("locusID_seed1", "bestSNP_lod_seed1", "locus_tag_seed1", 
    "locusID_seed2", "bestSNP_lod_seed2", "locus_tag_seed2", 
    "locusID_leaf2", "bestSNP_lod_leaf2", "locus_tag_leaf2",
    "met_rep1", "mz_rep1", "rt_rep1",
    "met_rep2", "mz_rep2", "rt_rep2")

for (i in 1:nrow(relat)) {
    print(paste(i, "/", nrow(relat)))
    
    ## relat[i, "seed1"] is name of metabolite in seed1
    inds_s1 <- which(seed1 == relat[i, "seed1"], arr.ind = TRUE)[, "row"] 
    locusID_inds_s1 <- seed1[inds_s1, "locusID"]
    locusID_inds_s1 <- unique(locusID_inds_s1)
    seed_inds1 <- seed1[seed1[, "locusID"] %in% locusID_inds_s1,]
    
    ## relat[i, "seed2"] is name of metabolite in seed2 corr to met in seed1
    inds_s2 <- which(seed2 == relat[i, "seed2"], arr.ind = TRUE)[, "row"] 
    locusID_inds_s2 <- seed2[inds_s2, "locusID"]
    locusID_inds_s2 <- unique(locusID_inds_s2)
    seed_inds2 <- seed2[seed2[, "locusID"] %in% locusID_inds_s2,]
    
    ## relat[i, "leaf2"] is name of metabolite in leaf2 corr to met in seed1
    inds_l2 <- which(leaf2 == relat[i, "leaf"], arr.ind = TRUE)[, "row"] 
    locusID_inds_l2 <- leaf2[inds_l2, "locusID"]
    locusID_inds_l2 <- unique(locusID_inds_l2)
    leaf_inds2 <- leaf2[leaf2[, "locusID"] %in% locusID_inds_l2,]
    
    ## write locus_tag in list: for each experiment in a list, in each 
    ## experiment list write a list for each locusID
    s1_l_lt <- list()
    s2_l_lt <- list()
    l2_l_lt <- list()
    
    ## write locus tags to a list 
    if (length(locusID_inds_s1)) {
        for (j in seq_along(locusID_inds_s1)) {
            s1_j <- seed_inds1[seed_inds1[, "locusID"] == locusID_inds_s1[j], ]
            s1_j_lt <- s1_j[, "locus_tag"]
            s1_l_lt[[j]] <- s1_j_lt
        }
    }
    
    if (length(locusID_inds_s2)) {
        for (j in seq_along(locusID_inds_s2)) {
            s2_j <- seed_inds2[seed_inds2[, "locusID"] == locusID_inds_s2[j], ]
            s2_j_lt <- s2_j[, "locus_tag"]
            s2_l_lt[[j]] <- s2_j_lt
        }
    }
    
    if (length(locusID_inds_l2)) {
        for (j in seq_along(locusID_inds_l2)) {
            l2_j <- leaf_inds2[leaf_inds2[, "locusID"] == locusID_inds_l2[j], ]
            l2_j_lt <- l2_j[, "locus_tag"]
            l2_l_lt[[j]] <- l2_j_lt
        }
    }
    
    ## get overlap and add additional information locusID, best_SNP_lod, locus_tag
    ol_l_i <- overlap_l(list(l1 = s1_l_lt, l2 = s2_l_lt, l3 = l2_l_lt))
    if (!is.matrix(ol_l_i)) {ol_l_i <- matrix(ol_l_i, ncol = 9)}
    
    ## get information on metabolites
    rel_rep1_i <- rel_rep1[paste("MLM.", rel_rep1[, 1], sep = "") == relat[i, 1], ]
    rel_rep2_i <- rel_rep2[paste("MLM.", rel_rep2[, 1], sep = "") == relat[i, 2], ]
    
    ## cbind to ol_l_i
    rep_t <- nrow(ol_l_i)
    res_i <- cbind(ol_l_i,  
        met_rep1 = rep(unique(rel_rep1_i[, "Name"]), rep_t),
        mz_rep1 = rep(unique(rel_rep1_i[, "mz"]), rep_t),
        rt_rep1 = rep(unique(rel_rep1_i[, "RT"]), rep_t),
        met_rep2 = rep(unique(rel_rep2_i[, "Name"]), rep_t),
        mz_rep2 = rep(unique(rel_rep2_i[, "mz"]), rep_t), 
        rt_rep2 = rep(unique(rel_rep2_i[, "RT"]), rep_t))

    ## rbind to res
    res <- rbind(res, res_i)
}


if (polarity == "neg") {
    write.table(res, file = "gwas_complete_met_all_neg.txt", sep = "\t", 
        dec = ".", quote = FALSE, row.names = FALSE)
} 
if (polarity == "pos") {
    write.table(res, file = "gwas_complete_met_all_pos.txt", sep = "\t", 
        dec = ".", quote = FALSE, row.names = FALSE)
}


## gwas_compete_met_all_neg/pos contains the best LOD per loci, not per metabolite,
## change to metabolite
options(stringsAsFactors = FALSE)

if (polarity == "neg") {
    setwd("~/GitHub/GWAS_arabidopsis_seed/")
    res <- read.table("gwas_complete_met_all_neg.txt", sep = "\t", dec = ".", header = TRUE) 
    setwd("~/GitHub/GWAS_arabidopsis_seed/results_build_loci")
    seed1 <- read.csv("gwas_locus_info_seed_1_neg.txt", header = TRUE, sep = "\t")
    seed2 <- read.csv("gwas_locus_info_seed_2_neg.txt", header = TRUE, sep = "\t")
    leaf2 <- read.csv("gwas_locus_info_leaf_neg.txt", header = TRUE, sep = "\t")
}
if (polarity == "pos") {
    setwd("~/GitHub/GWAS_arabidopsis_seed/")
    res <- read.table("gwas_complete_met_all_pos.txt", sep = "\t", dec = ".", header = TRUE) 
    setwd("~/GitHub/GWAS_arabidopsis_seed//results_build_loci")
    seed1 <- read.csv("gwas_locus_info_seed_1_pos.txt", header = TRUE, sep = "\t")
    seed2 <- read.csv("gwas_locus_info_seed_2_pos.txt", header = TRUE, sep = "\t")
    leaf2 <- read.csv("gwas_locus_info_leaf_pos.txt", header = TRUE, sep = "\t")
}

for (i in 1:nrow(res)) {
    locusID_seed1 <- res[i, "locusID_seed1"]
    locusID_seed2 <- res[i, "locusID_seed2"]
    locusID_leaf2 <- res[i, "locusID_leaf2"]
    met1_i <- paste("MLM.", res[i, "met_rep1"], sep = "")
    met2_i <- paste("MLM.", res[i, "met_rep2"], sep = "")
    
    if (!is.na(locusID_seed1)) {
        seed1_i <- seed1[seed1[, "locusID"] %in% locusID_seed1 & seed1[, "Peak.ID"] == met1_i, ]
        seed1_i_lod <- max(seed1_i[, "lod"])
        res[i, "bestSNP_lod_seed1"] <- seed1_i_lod
    }
    if (!is.na(locusID_seed2)) {
        seed2_i <- seed2[seed2[, "locusID"] %in% locusID_seed2 & seed2[, "Peak.ID"] == met2_i, ]
        seed2_i_lod <- max(seed2_i[, "lod"])
        res[i, "bestSNP_lod_seed2"] <- seed2_i_lod
    }
    if (!is.na(locusID_leaf2)) {
        leaf2_i <- leaf2[leaf2[, "locusID"] %in% locusID_leaf2 & leaf2[, "Peak.ID"] == met2_i, ]
        leaf2_i_lod <- max(leaf2_i[, "lod"])
        res[i, "bestSNP_lod_leaf2"] <- leaf2_i_lod
    }
}
if (polarity == "neg") {
    setwd("~/GitHub/GWAS_arabidopsis_seed/")
    write.table(res, file = "gwas_complete_met_all_trueLociLOD_neg.txt", 
        sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
} 
if (polarity == "pos") {
    setwd("~/GitHub/GWAS_arabidopsis_seed/")
    write.table(res, file = "gwas_complete_met_all_trueLociLOD_pos.txt", 
        sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
}