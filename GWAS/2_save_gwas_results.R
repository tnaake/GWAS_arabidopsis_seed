## R-script to create GWAS result file

p_value_thresh <- 0.001
path <- "."
work.dir <- "~/AG-Fernie/Thomas/Data/arabidopsis_seed/00_two_biological_replicates_results/"

##########################################

getResults <- function(path) {
    files <- dir(path, pattern = "GWAS.Results.csv", recursive = TRUE)
    results  <- list()
    for(f in files) {
        message(f)
        z <- regexpr("GAPIT\\.\\.?(.*)\\.GWAS.Results", f, perl = TRUE)
        stopifnot(z != -1)
        peak_id <- substring(f, attr(z, 'capture.start')[1],
                    attr(z, 'capture.start')[1] + attr(z, 'capture.length')[1] - 1)
        res <- read.csv(file.path(path, f))
        res$lod <- -log10(res$P.value)
        res$SNP <-  as.numeric( sub("m", "", res$SNP) )
        res$Peak.ID <- peak_id
        results[[peak_id]] <- res[res$P.value < p_value_thresh, ]
    }
    return(results)
}

## seed 1 positive
setwd("~/AG-Fernie/Thomas/Data/arabidopsis_seed/00_two_biological_replicates_results/seed_1_pos/")
gwas <- getResults(path)
save(gwas, file = file.path(work.dir, "./results_save_gwas_results/gwas_results_seed_1_pos.RData"))

## seed 1 negative
setwd("~/AG-Fernie/Thomas/Data/arabidopsis_seed/00_two_biological_replicates_results/seed_1_neg/")
gwas <- getResults(path)
save(gwas, file = file.path(work.dir, "./results_save_gwas_results/gwas_results_seed_1_neg.RData"))

## seed 2 positive
setwd("~/AG-Fernie/Thomas/Data/arabidopsis_seed/00_two_biological_replicates_results/seed_2_pos/")
gwas <- getResults(path)
save(gwas, file = file.path(work.dir, "./results_save_gwas_results/gwas_results_seed_2_pos.RData"))

## seed 2 negative
setwd("~/AG-Fernie/Thomas/Data/arabidopsis_seed/00_two_biological_replicates_results/seed_2_neg/")
gwas <- getResults(path)
save(gwas, file = file.path(work.dir, "./results_save_gwas_results/gwas_results_seed_2_neg.RData"))

## leaf positive
setwd("~/AG-Fernie/Thomas/Data/arabidopsis_seed/00_two_biological_replicates_results/leaf_pos/")
gwas <- getResults(path)
save(gwas, file = file.path(work.dir, "./results_save_gwas_results/gwas_results_leaf_pos.RData"))

## leaf negative 
setwd("~/AG-Fernie/Thomas/Data/arabidopsis_seed/00_two_biological_replicates_results/leaf_neg/")
gwas <- getResults(path)
save(gwas, file = file.path(work.dir, "./results_save_gwas_results/gwas_results_leaf_neg.RData"))

## explanation:
## the output results is a list. the number of "gwas" is the number of traits which is the number of metabolites here.
## for each object in the gwas stands for a metabolite, and there is another list for each metabolite containing 11 objects
## from 1-9 is coming from the result table, the 10th is LOD value
## the 11th is the Peak.ID for this metabolite
