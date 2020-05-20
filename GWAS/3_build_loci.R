## function to build loci
build_loci <- function(gwas = gwas, lod.thr = 5.3, pos.thr = 1e4) {
    z <- do.call('rbind', gwas)
    z <- z[order(z[,'Chromosome'], z[,'Position']), ]
    z <- z[z[,'lod'] > lod.thr, ]

    g <- rep(1, nrow(z))
    for(i in 2:nrow(z)) {
        if(z[i, 'Chromosome'] == z[i-1, 'Chromosome'] &
            abs(z[i, 'Position'] - z[i-1, 'Position'] ) < pos.thr) {
            g[i] <- g[i-1]
        } else {
            g[i] <- g[i-1] + 1
        }
    }
    print("g")
    rownames(z) <- NULL
    z <- data.frame(z)
    z$locusID <- g

    ## get genes
    genes <- lapply(unique(z$locusID), function(g) {
        w <- z[z$locusID == g, ]
        chr <- w$Chromosome[1]
        x  <- min(w$Pos)
        y  <- max(w$Pos)
        j <- which(chr == tair$chrom & x - tair$pos_f < pos.thr & tair$pos_i - y < pos.thr)
        tair$locus_tag[j]
    })
    print("genes")

    ## merge common genes and reassign locus ID
    merge  <- list(1)
    k <- 1
    for(i in 2:length(genes)) {
        if(any(genes[[i]] %in% unlist(genes[ merge[[k]] ]))) {
            merge[[k]] <- c( merge[[k]] , i )
        } else {
            k <- k + 1
            merge[[k]] <- i
        }
    }
    print("k")

    tmp <- list()
    tmpID <- rep(0, nrow(z))
    for(i in 1:length(merge)) {
        tmp[[i]] <- unique(unlist(genes[merge[[i]]]))
        tmpID[ z$locusID %in% merge[[i]] ] <- i
    }

    genes <- tmp
    stopifnot(length(unlist(genes)) == length(unique(unlist(genes))))
    z$locusID <- tmpID
    
    z <- z[order(z$locusID, -z$lod), ]
    loci_nfo <- z
    
    ngenes <- sapply(genes, length)
    
    tmp <- unlist(genes)
    tmp2 <- tair[match(tmp, tair$locus_tag), c('chrom', 'pos_i', 'pos_f')]
    rownames(tmp2) <- tmp
    
    # compute distances between genes
    D <- matrix(NA, nrow(tmp2), nrow(z))
    rownames(D) <- tmp

    for(i in 1:nrow(tmp2)) {
        for(j in 1:nrow(z)) {
            if(z$Chromosome[j] != tmp2$chrom[i])
                next
            if(z$Position[j] >= tmp2$pos_f[i])
                D[i, j] <- tmp2$pos_f[i] - z$Position[j]
            else if(z$Position[j] <= tmp2$pos_i[i])
                D[i, j] <- tmp2$pos_i[i] - z$Position[j]
            else
                D[i, j] <- 0
        }
    }
    print("D")

    nearest <- list()
    best <- list()
    for(i in 1:length(genes)) {
        j <- which(z$locusID == i)
        best[[i]] <- rep(j[1], length(genes[[i]]))
        d <- D[genes[[i]], j, drop = F]
        nearest[[i]] <- j[apply(d, 1, function(w) which.min(abs(w)))]
    }
    
    nearest <- unlist(nearest)
    best <- unlist(best)
    
    gene_nfo <- data.frame(
        locusID = rep(1:length(genes), times = ngenes),
        no_SNP = rep(table(z$locusID), times = ngenes),
        best_SNP = z$SNP[best],
        best_SNP_lod = z$lod[best],
        best_SNP_dist = mapply(function(i, j) D[i, j], 1:length(best), best),
        nearest_SNP = z$SNP[nearest],
        nearest_SNP_lod = z$lod[nearest],
        nearest_SNP_dist = mapply(function(i, j) D[i, j], 1:length(nearest), nearest),
        best_eq_nearest = best == nearest
    )
    print("gene_nfo")
    tr = tair[match(tmp, tair$locus_tag), ]
    gene_nfo <- cbind(gene_nfo, tr)
    
    # this section add the metabolite names to the gene_info data.
    z <- split(loci_nfo$Peak.ID, loci_nfo$locusID)
    z <- sapply(z, unique, simplify = FALSE)
    m <- max(sapply(z, length))
    
    if(m == 1) {
        z <- cbind(Peak.ID = unlist(z))
    } else {
        z <- t(sapply(z, function(w) if(length(w) == m) w else c(w, rep(".", m - length(w)))))
        colnames(z) <- sprintf("Peak.ID_%s", 1:m)
    }
    z <- z[gene_nfo$locusID,, drop = FALSE]
    gene_nfo <- cbind(gene_nfo, z)

    return(list(loci_nfo, gene_nfo))
}

##load("~/winhome/Documents/01_GWAS/01_Data/00_two_biological_replicates/TAIR9.RData")
load("/home/naake/01_GWAS/TAIR9.RData")
## LOD threshold
lod.thr <- 5.3
## window in base pairs (plus or minus)
pos.thr <- 1e4 ## 2e4

## set working directory
##setwd("~/AG-Fernie/Thomas/Data/arabidopsis_seed/00_two_biological_replicates_results")
setwd("/home/naake/01_GWAS/")
## leaf negative
# load("./results_save_gwas_results/gwas_results_leaf_neg.RData")
# leaf_neg <- build_loci(gwas, lod.thr = lod.thr, pos.thr = pos.thr)
# write.table(leaf_neg[[1]], file = './results_build_loci/gwas_locus_info_leaf_neg.txt', sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(leaf_neg[[2]], file = './results_build_loci/gwas_gene_info_leaf_neg.txt', sep = "\t", quote = FALSE, row.names = FALSE)

## leaf positive
# load("./results_save_gwas_results/gwas_results_leaf_pos.RData")
# leaf_pos <- build_loci(gwas, lod.thr = lod.thr, pos.thr = pos.thr)
# write.table(leaf_pos[[1]], file = './results_build_loci/gwas_locus_info_leaf_pos.txt', sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(leaf_pos[[2]], file = './results_build_loci/gwas_gene_info_leaf_pos.txt', sep = "\t", quote = FALSE, row.names = FALSE)

## seed 1 negative
# load("./results_save_gwas_results/gwas_results_seed_1_neg.RData")
# seed_1_neg <- build_loci(gwas, lod.thr = lod.thr, pos.thr = pos.thr)
# write.table(seed_1_neg[[1]], file = './results_build_loci/gwas_locus_info_seed_1_neg.txt', sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(seed_1_neg[[2]], file = './results_build_loci/gwas_gene_info_seed_1_neg.txt', sep = "\t", quote = FALSE, row.names = FALSE)

## seed 1 positive
# load("./results_save_gwas_results/gwas_results_seed_1_pos.RData")
# seed_1_pos <- build_loci(gwas, lod.thr = lod.thr, pos.thr = pos.thr)
# write.table(seed_1_pos[[1]], file = './results_build_loci/gwas_locus_info_seed_1_pos.txt', sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(seed_1_pos[[2]], file = './results_build_loci/gwas_gene_info_seed_1_pos.txt', sep = "\t", quote = FALSE, row.names = FALSE)

## seed 2 negative
# load("./results_save_gwas_results/gwas_results_seed_2_neg.RData")
# seed_2_neg <- build_loci(gwas, lod.thr = lod.thr, pos.thr = pos.thr)
# write.table(seed_2_neg[[1]], file = './results_build_loci/gwas_locus_info_seed_2_neg.txt', sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(seed_2_neg[[2]], file = './results_build_loci/gwas_gene_info_seed_2_neg.txt', sep = "\t", quote = FALSE, row.names = FALSE)

## seed 2 positive
load("./results_save_gwas_results/gwas_results_seed_2_pos.RData")
seed_2_pos <- build_loci(gwas, lod.thr = lod.thr, pos.thr = pos.thr)
write.table(seed_2_pos[[1]], file = './results_build_loci/gwas_locus_info_seed_2_pos.txt', sep = "\t", quote = FALSE, row.names = FALSE)
write.table(seed_2_pos[[2]], file = './results_build_loci/gwas_gene_info_seed_2_pos.txt', sep = "\t", quote = FALSE, row.names = FALSE)