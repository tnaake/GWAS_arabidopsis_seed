## load SNP data
setwd("~/Documents/01_GWAS/01_Data/GWAS_script/LD/")
load("../LD/Data/TAIR9.RData")

## load package ape
## ape is used for calculating the distance between different haplotypes
library(ape) 

## input the association to check
setwd("W:/Thomas/Data/arabidopsis_seed/00_two_biological_replicates_results/results_build_loci")
options(stringsAsFactors=FALSE)
PC_input_1 <- read.table("gwas_gene_info_seed_1_neg.txt", header = TRUE, 
    sep = "\t", quote = "\"", dec = ".", comment.char = "")
PC_input_2 <- read.table("gwas_gene_info_seed_2_neg.txt", header = TRUE, 
    sep = "\t", quote = "\"", dec = ".", comment.char = "")

## genes to test:
## enter here a character vector of genes to test
goi <- c("AT5G17040", "AT5G17050")
# goi <- c(
#     "AT4G03050", ## AOP3
#     "AT4G03060", ## AOP2
#     "AT4G03063", ## AOP of AT4G03070; AOP1 pseudogene
#     "AT4G03070", ## AOP1
#     "AT5G23010", ## MAM1
#     "AT5G23030") ## MAM3
    

## metabolites to test: 
## enter here the matched mass feature pairs that correspond to the metabolites
## in the form of a two column matrix
moi <- matrix(c("Cluster_28337", "Cluster_057215",
    "Cluster_35169", "Cluster_075199",
    "Cluster_36050", "Cluster_077251"),
    ncol = 2, byrow = TRUE)
# moi <- matrix(c("Cluster_27139", "Cluster_053358", ## 4-methylsulfinylbutyl GS
#   "Cluster_27723", "Cluster_055415", ## 5-methylsulfinylpentyl GS
#   "Cluster_28387", "Cluster_057362", ## 6-methylsulfinylhexyl GS
#   "Cluster_29158", "Cluster_059324", ## 7-methylsulfinylheptyl GS
#   "Cluster_29750", "Cluster_061218", ## 8-methylsulfinyloctyl GS
#   "Cluster_24564", "Cluster_044614", ## 3-hydroxypropyl GS
#   "Cluster_25248", "Cluster_046712", ## 4-hydroxybutyl GS
#   "Cluster_25869", "Cluster_049033", ## 3-methylthiopropyl GS
#   "Cluster_26400", "Cluster_051010", ## 4-methylthiobutyl GS
#   "Cluster_27077", "Cluster_053100", ## 5-methylthiopentyl GS
#   "Cluster_28299", "Cluster_057117", ## 7-methylthipheptyl GS
#   "Cluster_29064", "Cluster_059172", ## 8-methylthiooctyl GS
#   "Cluster_29119", "Cluster_059172", ## 1-methoxy-3-indolylmethyl GS
#   "Cluster_29269", "Cluster_059592", ## 3-benzoyloxypropyl GS
#   "Cluster_29899", "Cluster_061425", ## 4-benzoyloxybutyl GS
#   "Cluster_30574", "Cluster_063305", ## 5-benzoyloxypentyl GS
#   "Cluster_31267", "Cluster_065294", ## 6-benzoyloxybenzyl GS
#   "Cluster_34455", "Cluster_073649", ## 3-sinapoyloxypropyl GS
#   "Cluster_02069", "Cluster_005378"),## 4-sinapoyloxybutyl GS
#   ncol = 2, byrow = TRUE)

## truncate the input
PC_input_1 <- cbind(PC_input_1[, "locusID"], PC_input_1[, "best_SNP_lod"], 
                PC_input_1[, "locus_tag"], PC_input_1[, 27:ncol(PC_input_1)])
PC_input_2 <- cbind(PC_input_2[, "locusID"], PC_input_2[, "best_SNP_lod"], 
                PC_input_2[, "locus_tag"], PC_input_2[, 27:ncol(PC_input_2)])

colnames(PC_input_1)[1:3] <- c("locusID", "LOD", "Gene")
colnames(PC_input_2)[1:3] <- c("locusID", "LOD", "Gene")

PC_input_1 <- PC_input_1[PC_input_1[, "Gene"] %in% goi, ]
PC_input_2 <- PC_input_2[PC_input_2[, "Gene"] %in% goi, ]

PC_input_1 <- data.frame(
    locusID = rep(PC_input_1[, "locusID"], length(moi[, 1])), 
    LOD = rep(PC_input_1[, "LOD"], length(moi[, 1])), 
    expand.grid(PC_input_1[, "Gene"], moi[, 1]))
PC_input_2 <- data.frame(
    locusID = rep(PC_input_2[, "locusID"], length(moi[, 2])), 
    LOD = rep(PC_input_2[, "LOD"], length(moi[, 2])), 
    expand.grid(PC_input_2[, "Gene"], moi[, 2]))

colnames(PC_input_1)[3:4] <- c("Gene", "metabolite")
colnames(PC_input_2)[3:4] <- c("Gene", "metabolite")

## load metabolite data
Data_1 <- read.table("~/Documents/01_GWAS/01_Data/00_two_biological_replicates/DMN_ecotypes_seed1_negative.txt", 
    header = TRUE, sep = "\t",quote = "\"",dec = ".",comment.char = "")
Data_2 <- read.table("~/Documents/01_GWAS/01_Data/00_two_biological_replicates/DMN_ecotypes_seed2_negative.txt", 
    header = TRUE, sep = "\t",quote = "\"",dec = ".",comment.char = "")

## remove m/z and retention time
Data_1 <- Data_1[-c(1:2), ]
Data_2 <- Data_2[-c(1:2), ]

## remove duplicated accessions
Data_1 <- Data_1[!duplicated(Data_1[, 1]), ]
Data_2 <- Data_2[!duplicated(Data_2[, 1]), ]

## actual accessions (remove the accessions which are absent in LC-MS dataset)
Taxa_1 <- as.character(Data_1[, 1])
Taxa_2 <- as.character(Data_2[, 1])
Taxa <- intersect(Taxa_1, Taxa_2)

## load SNP info (contains SNP from all chromosomes):
load("~/Documents/01_GWAS/01_Data/GWAS_script/LD/Data/SNP_chr2.RData")
SNP_chr2_part1 <- SNP_chr2[, 1:4]
SNP_chr2_part2 <- SNP_chr2[ , 12:ncol(SNP_chr2)]

SNP_chr2_part2_new <- SNP_chr2_part2[, colnames(SNP_chr2_part2) %in% Taxa]
SNP_chr2_part1_new <- SNP_chr2_part1[, -2]

SNP_chr2 <- cbind(SNP_chr2_part1_new, SNP_chr2_part2_new) ## dim=199455*number of accessions

## function for clustering SNPs in Gene
Cluster_Part <- function(Gene) {
    a <- which(tair[[6]] == Gene)
    Chr_Nr <- tair[[1]][a]
    pos_start <- tair[[2]][a]
    pos_end <- tair[[3]][a]
    rowNr <- 1:nrow(SNP_chr2)
    a <- which(SNP_chr2[, 2] == Chr_Nr)
    SNP_NO <- rowNr[(SNP_chr2[, 3] > (pos_start)) & (SNP_chr2[, 3] < (pos_end))]
    ## check the SNPs picked up: SNP_chr2[(SNP_NO[SNP_NO %in% a]),1:4]
    SNP_NO <- SNP_NO[SNP_NO %in% a]
    
    if (length(SNP_NO) < 2) {
        a <- which(tair[[6]] == Gene)
        Chr_Nr <- tair[[1]][a]
        pos_start <- tair[[2]][a]
        pos_end <- tair[[3]][a]
        dis1 <- pos_start-tair[[3]][a - 1]
        dis2 <- tair[[2]][a + 1] - pos_end
        if (dis1 < 0) {print("Attention!wrong calculation!")}
        if (dis2 < 0) {print("Attention!wrong calculation!")}
        rowNr <- 1:nrow(SNP_chr2)
        a <- which(SNP_chr2[, 2] == Chr_Nr)
        SNP_NO <- rowNr[ (SNP_chr2[, 3] > (pos_start - dis1 + 1)) & 
                                        (SNP_chr2[, 3] < (pos_end + dis2 - 1)) ]
        ## check the SNPs picked up: SNP_chr2[(SNP_NO[SNP_NO %in% a]), 1:4]
        SNP_NO <- SNP_NO[SNP_NO %in% a]
    }

    if (length(SNP_NO) < 2) {
        a <- which(tair[[6]] == Gene)
        Chr_Nr <- tair[[1]][a]
        pos_start <- tair[[2]][a]
        pos_end <- tair[[3]][a]
        dis <- 1000
        rowNr <- 1:nrow(SNP_chr2)
        a <- which(SNP_chr2[, 2] == Chr_Nr)
        SNP_NO <- rowNr[ (SNP_chr2[, 3] > (pos_start - dis)) & 
                                        (SNP_chr2[, 3] < (pos_end + dis)) ]
        
        ## check the SNPs picked up: SNP_chr2[(SNP_NO[SNP_NO %in% a]),1:4]
        SNP_NO <- SNP_NO[SNP_NO %in% a]
    }

    if (length(SNP_NO) < 2) {
        a <- which(tair[[6]] == Gene)
        Chr_Nr <- tair[[1]][a]
        pos_start <- tair[[2]][a]
        pos_end <- tair[[3]][a]
        dis <- 1500
        rowNr <- 1:nrow(SNP_chr2)
        a <- which(SNP_chr2[, 2] == Chr_Nr)
        SNP_NO <- rowNr[ (SNP_chr2[, 3] > (pos_start-dis)) & 
                                        (SNP_chr2[, 3] < (pos_end + dis)) ]
        ## check the SNPs picked up: SNP_chr2[(SNP_NO[SNP_NO %in% a]),1:4]
        SNP_NO <- SNP_NO[SNP_NO %in% a]
    }

    if (length(SNP_NO) < 2) {
        a <- which(tair[[6]] == Gene)
        Chr_Nr <- tair[[1]][a]
        pos_start <- tair[[2]][a]
        pos_end <- tair[[3]][a]
        dis <- 2000
        rowNr <- 1:nrow(SNP_chr2)
        a <- which(SNP_chr2[, 2] == Chr_Nr)
        SNP_NO <- rowNr[(SNP_chr2[, 3]>(pos_start - dis)) & 
                                        (SNP_chr2[, 3] < (pos_end + dis))]
        ## check the SNPs picked up: SNP_chr2[(SNP_NO[SNP_NO %in% a]),1:4]
        SNP_NO <- SNP_NO[SNP_NO %in% a]
    }

    if (length(SNP_NO) < 2) {
        a <- which(tair[[6]] == Gene)
        Chr_Nr <- tair[[1]][a]
        pos_start <- tair[[2]][a]
        pos_end <- tair[[3]][a]
        dis <- 4000
        rowNr <- 1:nrow(SNP_chr2)
        a <- which(SNP_chr2[,2] == Chr_Nr)
        SNP_NO <- rowNr[(SNP_chr2[, 3]>(pos_start - dis)) & 
                                        (SNP_chr2[, 3] < (pos_end + dis))]
        ## check the SNPs picked up: SNP_chr2[(SNP_NO[SNP_NO %in% a]),1:4]
        SNP_NO <- SNP_NO[SNP_NO %in% a]
    }

    haplotype_all <- NULL
    for (i in 4:ncol(SNP_chr2)) {
        b <- ""
        for (j in 1:length(SNP_NO)) {
            b <- paste(b, as.character(SNP_chr2[SNP_NO[j], i]), sep = "")
        }
        haplotype_all <- c(haplotype_all, b)
    }

    haplotype <- unique(haplotype_all)
    ## SNP name and SNP info cross accessions
    haplo <- SNP_chr2[SNP_NO, -(2:3)]  
    test_haplo <- haplo[, -1]
    rownames(test_haplo) <- as.character(haplo[, 1])
    haplo_forDist <- t(test_haplo)

    ## calculate distance between haplotypes
    x1 <- dist.gene(haplo_forDist, method = "pairwise", 
        pairwise.deletion = FALSE, variance = FALSE)
    hclust1 <- hclust(x1, "ward.D")
    memb <- cutree(hclust1, h = 0.00001)
    ## simplify the tree
    haplo_forDist_simple <- NULL
    rownames_haploS <- NULL
    for (i in 1:length(unique(memb))) {
        a <- (which(memb == unique(memb)[i]))[1]
        rownames_haploS <- c(rownames_haploS, rownames(haplo_forDist)[a])
        haplo_forDist_simple <- rbind(haplo_forDist_simple, haplo_forDist[a, ])
    }
    rownames(haplo_forDist_simple) <- rownames_haploS
    x2 <- dist.gene(haplo_forDist_simple, method = "pairwise", 
        pairwise.deletion = FALSE, variance = FALSE)
    hclust2 <- hclust(x2, "ward.D")
    hclust2$labels <- as.character(table(memb))
    memb1 <- cutree(hclust2, h = 0.00001)
    output_cluster <- list()
    output_cluster[[1]] <- hclust2
    output_cluster[[2]] <- memb
    output_cluster[[3]] <- hclust1
    return(output_cluster)
}

### box-plot
# New function for boxplot
Boxplot_Part <- function(memb, hclust2, Data) {
    haplo_Type <- unique(memb)
    hclust2$order
    haplo_name <- NULL
    for (i in 1:length(hclust2$order)) {
        a <- paste("H", hclust2$order[i], sep = "")
        aa <- rep(a, length(which(memb == hclust2$order[i])))
        cc <- 1:length(which(memb == hclust2$order[i]))
        b <- cbind(which(memb == hclust2$order[i]), aa, cc)
        haplo_name <- rbind(haplo_name, b)
    }
    haplo_name_new <- haplo_name[, -1]
    colnames(haplo_name_new) <- c("sample", "replica")

    ## prepare sample list file
    samp <- haplo_name_new  

    Data_1 <- Data[, 2:ncol(Data)]
    rownames(Data_1) <- as.character(Data[, 1])
    Data_2 <- t(Data_1)
    Data_range <- NULL
    for (i in 1:nrow(samp)) {
        a <- which(colnames(Data_2) == rownames(samp)[i])
        Data_range <- cbind(Data_range, Data_2[, a])
    }
    colnames(Data_range) <- rownames(samp)
    metab <- Data_range
    groups <- factor(as.factor(haplo_name_new[, 1]), 
                                levels = unique(as.factor(haplo_name_new[, 1])))
    group.length <- length(levels(groups))
    output_boxplot <- list()
    output_boxplot[[1]] <- metab
    output_boxplot[[2]] <- groups
    output_boxplot[[3]] <- group.length
    return(output_boxplot)
}

## Function to get a list of all group-pairwise combinations
getPairs <- function(g) {
    z <- levels(g)
    out <- character(length(z) * (length(z) - 1) / 2)
    k <- 1
    for(i in 1:(length(z) - 1))
        for(j in (i+1):length(z)) {
            out[k] <- paste(z[j], z[i], sep = "-")
            k <- k + 1
        }
    return(out)
}

## start here with the actual analysis
## run the following for replicate 1 and replicate 2
setwd("~/Documents/01_GWAS/01_Data/GWAS_script/haplotype/")
rep <- "rep_1"

if (rep == "rep_1") {
    PC_input <- PC_input_1
    Data <- Data_1
}
    
if (rep == "rep_2") {
    PC_input <- PC_input_2
    Data <- Data_2
}

for (i in 1:nrow(PC_input)) {
    print(i)
    MetID <- as.character(PC_input[i, "metabolite"])
    Gene <- as.character(PC_input[i, "Gene"])
    OUTPUT_clust <- Cluster_Part(Gene)
    hclust2 <- OUTPUT_clust[[1]]
    OUTPUT_boxplot <- Boxplot_Part(OUTPUT_clust[[2]], hclust2, Data)
    metab <- OUTPUT_boxplot[[1]]
    groups <- OUTPUT_boxplot[[2]]
    
    ## P-values
    anova.2 <- function(x, y) anova(lm(x ~ y))$Pr[1]
    p.values <- anova.2(metab[MetID, ], groups)
    p.values.adj <- p.adjust(p.values, "BH", n = nrow(metab))
    
    ## Tukey test (test on differences between the means of the levels of a factor)
    tukey.2 <- function(x, y) TukeyHSD(aov(x ~ y))
    tuk <- tukey.2(metab[MetID, ], groups)

    gp <- getPairs(groups)
    ## extract the TukeyHSD p-values
    tuk.p.values <- tuk$y[, "p adj"]
    tup.p.values <- tuk.p.values[gp]
    names(tuk.p.values) <- gp

    pdf(paste(i,"_", MetID, "_AND_", Gene,".pdf", sep=""), 12, 8)
    par(fig=c(0, 1, 0, 0.7), new = TRUE)
    boxplot(metab[MetID, ] ~ groups, col = rainbow(OUTPUT_boxplot[[3]]), las = 1,
	    main = NULL, ylab = "Intensity",
		sub = sprintf("ANOVA: p-value = %.2e (FDR corrected)", p.values.adj), 
		frame = FALSE)
    par(fig = c(0.025, 0.975, 0.5, 1), new = TRUE)
    ## create plot
    plot(hclust2, hang = -1, sub = "", xlab = "",
         main = paste("Gene:", Gene," | Metabolite:", MetID, sep = ""))
    dev.off()
}
