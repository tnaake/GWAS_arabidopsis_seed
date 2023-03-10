### Calculate broad-sense heritability:
options(stringsAsFactors = FALSE)
setwd("~/GitHub/GWAS_arabidopsis_seed/")
ion_mode <- "neg"

## for positive mode
if (ion_mode == "pos") {
  
    ## load data with matching information
    rel_rep1 <- read.csv("rep1_positive_match.csv", sep="\t")
    rel_rep2 <- read.csv("rep2_positive_match.csv", sep="\t")
    
    ## load metabolite data
    met_rep1 <- read.table("DMN_ecotypes_seed1_positive.txt", header = TRUE)
    met_rep2 <- read.table("DMN_ecotypes_seed2_positive.txt", header = TRUE)
}

if (ion_mode == "neg") {
    
    ## load data with matching information
    rel_rep1 <- read.csv("rep1_negative_match.csv", sep="\t")
    rel_rep2 <- read.csv("rep2_negative_match.csv", sep="\t")
    
    ## load metabolite data
    met_rep1 <- read.table("DMN_ecotypes_seed1_negative.txt", header = TRUE)
    met_rep2 <- read.table("DMN_ecotypes_seed2_negative.txt", header = TRUE)
}

## filter based on rt_dev and mz_dev
inds_keep <- logical(nrow(rel_rep1))
inds_keep[rel_rep1[, "rt_dev"] <= 0.075 & abs(rel_rep1[, "mz_dev"]) <= 0.0075] <- TRUE
rel_rep1 <- rel_rep1[inds_keep, ]
rel_rep2 <- rel_rep2[inds_keep, ]

## remove duplicated taxa
met_rep1 <- met_rep1[!duplicated(met_rep1[, "Taxa"]), ]
met_rep2 <- met_rep2[!duplicated(met_rep2[, "Taxa"]), ]
rownames(met_rep1) <- met_rep1[, "Taxa"]
rownames(met_rep2) <- met_rep2[, "Taxa"]

## remove RT and mz
met_rep1 <- met_rep1[!rownames(met_rep1) %in% c("RT", "mz"), ]
met_rep2 <- met_rep2[!rownames(met_rep2) %in% c("RT", "mz"), ]
## remove Taxa
met_rep1 <- met_rep1[, -which(colnames(met_rep1) == "Taxa")]
met_rep2 <- met_rep2[, -which(colnames(met_rep2) == "Taxa")]

## remove all column pairs that have only zero
dim(rel_rep1)
dim(rel_rep2)
inds_rem1 <- which(apply(met_rep1, 2, function(x) all(x == 0)))
inds_rem2 <- which(apply(met_rep2, 2, function(x) all(x == 0)))
inds_rem1 <- which(rel_rep1[, "Name"] %in% names(inds_rem1)) ## indices in rel_rep1 to remove
inds_rem2 <- which(rel_rep2[, "Name"] %in% names(inds_rem2)) ## indices in rel_rep2 to remove
inds_rem <- unique(c(inds_rem1, inds_rem2))
all(rel_rep1[, "mapping_rep1_rep2"] == rel_rep2[, "mapping_rep1_rep2"])
rel_rep1 <- rel_rep1[-inds_rem, ]
rel_rep2 <- rel_rep2[-inds_rem, ]
dim(rel_rep1)
dim(rel_rep2)

## only retain those metabolic features that are in rel_rep1 and rep2
all(rel_rep1[, "mapping_rep1_rep2"] == rel_rep2[, "mapping_rep1_rep2"])
dim(met_rep1)
dim(met_rep2)
met_rep1 <- met_rep1[, rel_rep1[, "Name"]]
met_rep2 <- met_rep2[, rel_rep2[, "Name"]]
dim(met_rep1)
dim(met_rep2)
which(apply(met_rep1, 2, function(x) all(x == 0)))
which(apply(met_rep2, 2, function(x) all(x == 0)))

## match rownames
met_rep1 <- met_rep1[rownames(met_rep1) %in% rownames(met_rep2), ]
met_rep2 <- met_rep2[rownames(met_rep2) %in% rownames(met_rep1), ]
inds <- match(rownames(met_rep1), rownames(met_rep2))
all(rownames(met_rep1) == rownames(met_rep2)[inds])
met_rep2_new <- met_rep2[inds, ]

## z-scale
met_rep1 <- apply(met_rep1, 2, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
met_rep2_new <- apply(met_rep2_new, 2, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

colnames(met_rep2_new) <- colnames(met_rep1)
H2_input <- rbind(met_rep1, met_rep2_new)

## filter based on correlation > 0.1
## calculate correlation values per matches mass feature pairs
cor_values <- numeric()
for (i in 1:ncol(H2_input)) {
    if (ion_mode == "pos") {
        cor_values <- c(cor_values, cor(H2_input[1:295, i], H2_input[296:nrow(H2_input), i]))  
    }
    if (ion_mode == "neg") {
        cor_values <- c(cor_values, cor(H2_input[1:296, i], H2_input[297:nrow(H2_input), i])) 
    }
}
inds_keep <- cor_values > 0.1
met_rep1 <- met_rep1[, inds_keep]
met_rep2_new <- met_rep2_new[, inds_keep]
H2_input <- rbind(met_rep1, met_rep2_new)

## create relation between the two replicates
n <- nrow(met_rep1)
H2_head <- data.frame(LINE = rep(1:n, 2), LOC = c(rep(1, n), rep(2, n)))
H2_input_new <- cbind(H2_head, H2_input)
  
## assign id and x
id <- as.character(H2_input_new[, 1]) 
x <- H2_input_new

## calculate broad-sense heritability, H2
library(Matrix)
library(lattice)
library(lme4)

line.blup <- list()
heritability <- list()

for (i in 3:ncol(x)) {
    print(paste(i, "/", ncol(x)))
    varcomp <- tryCatch(lmer(x[, i] ~ (1 | LINE) + (1 | LOC), data = x), error = function(e) return(NULL))
    if (!is.null(varcomp)) {
        var.trans = lme4::VarCorr(varcomp)
        var = data.frame(Groups = c('LINE', 'LOC', 'Residual'), 
            Variance = c(as.numeric(var.trans$LINE),
                as.numeric(var.trans$LOC), attr(var.trans,'sc')^2), 
                check.names = FALSE)
        
        ## residual standard deviation is stored as attribute "sc" 
        gvar <- as.numeric(as.character(var$Variance))[var$Groups %in% 'LINE']
        evar <- as.numeric(as.character(var$Variance))[var$Groups %in% 'Residual']
        heritability[[i-2]] <- matrix(gvar / (gvar + evar / length(unique(x$LOC))), 1, 1, dimnames = list(colnames(x)[i], 'H2'))
        f <- fixef(varcomp)
        r <- ranef(varcomp)$LINE
        blup <- f + r  
        line.blup[[i-2]] <- blup[match(id, rownames(blup)), ]
        
    } else {
        print("Error!")
        heritability[[i-2]] <- matrix(0, 1, 1, dimnames = list(colnames(x)[i], 'H2'))
        line.blup[[i-2]] <- rep(NaN, nrow(x))
    }
}

heritability <- do.call("rbind", heritability)
line.blup <- cbind(id, do.call("cbind", line.blup))
colnames(line.blup) <- c("line", names(x)[-c(1:2)])

## write to files
setwd("~/GitHub/GWAS_arabidopsis_seed/heritability")

if (ion_mode == "pos") {
    write.table(heritability, "BroadSenseHeritability3e_pos_seed.csv", sep = ",", quote = F)
    write.csv(line.blup, "BlupEachlines3e_pos_seed.csv", row.names = FALSE)
    write.table(H2_input_new, "H2_input_pos.csv", sep = ",", row.names = FALSE)
}

if (ion_mode == "neg") {
    write.table(heritability, "BroadSenseHeritability3e_neg_seed.csv", sep = ",", quote = F)
    write.csv(line.blup, "BlupEachlines3e_neg_seed.csv", row.names = FALSE)
    write.table(H2_input_new, "H2_input_neg.csv", sep = ",", row.names = FALSE)
}

## create plot
library(ggplot2)

## load the data sets and results
options(stringsAsFactors = FALSE)
heritability_pos <- read.table("BroadSenseHeritability3e_pos_seed.csv", sep = ",", header = TRUE)
heritability_neg <- read.table("BroadSenseHeritability3e_neg_seed.csv", sep = ",", header = TRUE)
h2_input_pos <- read.table("H2_input_pos.csv", sep = ",", header = TRUE)
h2_input_neg <- read.table("H2_input_neg.csv", sep = ",", header = TRUE)

## get for mapped ones
mapped_pos <- read.table("../gwas_complete_met_all_pos.txt", sep = "\t", header = TRUE)
mapped_neg <- read.table("../gwas_complete_met_all_neg.txt", sep = "\t", header = TRUE)    

## mapped for both replicates
mapped_uniq_pos <- unique(mapped_pos[mapped_pos$locus_tag_seed1 != "" & 
    mapped_pos$locus_tag_seed2 != "", "met_rep1"])
inds_map_pos <- mapped_uniq_pos[mapped_uniq_pos %in% rownames(heritability_pos)]
hist(heritability_pos[, 1])
hist(heritability_pos[inds_map_pos, 1])

mapped_uniq_neg <- unique(mapped_neg[mapped_neg$locus_tag_seed1 != "" & 
    mapped_neg$locus_tag_seed2 != "", "met_rep1"])
inds_map_neg <- mapped_uniq_neg[mapped_uniq_neg %in% rownames(heritability_neg)]
hist(heritability_neg[,1])
hist(heritability_neg[inds_map_neg,])

## plot heritability
df <- data.frame(
    values = c(heritability_pos[, 1], heritability_pos[inds_map_pos,],
        heritability_neg[, 1], heritability_neg[inds_map_neg,]),
    ion_mode = c(rep("pos", length(c(heritability_pos[, 1], inds_map_pos))),
        rep("neg", length(c(heritability_neg[, 1], inds_map_neg)))),
    type = c(rep("all", length(heritability_pos[, 1])), rep("mapped", length(inds_map_pos)),
        rep("all", length(heritability_neg[, 1])), rep("mapped", length(inds_map_neg))))

dodge <- position_dodge(width = 0.6)
g <- ggplot(df, aes(y = values, x = type, fill = type)) + 
    geom_violin(width = 0.7, scale = "count", position = dodge, alpha = 0.3) + 
    geom_boxplot(width = 0.1, color = "black", alpha = 0.75, position = dodge) + facet_grid(~ion_mode) +
    theme_bw() + ylab("heritability") 

setwd("~/GitHub/GWAS_arabidopsis_seed/heritability")
ggsave(g, filename = "plot_heritability_mass_features_violin.pdf", device = "pdf")

## plot correlation values between matched metabolite pairs
cor_pos <- cor_neg <- numeric()
for (i in 1:ncol(h2_input_pos)) {
    cor_pos <- c(cor_pos, cov(h2_input_pos[1:295, i], h2_input_pos[296:nrow(h2_input_pos), i]))
}
for (i in 1:ncol(h2_input_neg)) {
    cor_neg <- c(cor_neg, cov(h2_input_neg[1:296, i], h2_input_neg[297:nrow(h2_input_neg), i]))
}
cor_pos_map <- cor_pos[match(inds_map_pos, colnames(h2_input_pos))]
cor_neg_map <- cor_neg[match(inds_map_neg, colnames(h2_input_neg))]
cor_pos <- cor_pos[-which.max(cor_pos)]
cor_neg <- cor_neg[-which.max(cor_neg)]

df <- data.frame(
    values = c(cor_pos, cor_pos_map, cor_neg, cor_neg_map),
    ion_mode = c(rep("pos", length(c(cor_pos, cor_pos_map))),
        rep("neg", length(c(cor_neg, cor_neg_map)))),
    type = c(rep("all", length(cor_pos)), rep("mapped", length(cor_pos_map)),
        rep("all", length(cor_neg)), rep("mapped", length(cor_neg_map))))

dodge <- position_dodge(width = 0.6)
g <- ggplot(df, aes(y = values, x = type, fill = type)) + 
  geom_violin(width = 0.7, scale = "count", position = dodge, alpha = 0.3) + 
  geom_boxplot(width = 0.1, color = "black", alpha = 0.75, position = dodge) + facet_grid(~ion_mode) +
  theme_bw() + ylab("Pearson correlation") 

setwd("~/GitHub/GWAS_arabidopsis_seed/heritability")
ggsave(g, filename = "plot_correlation_mass_features_violin.pdf", device = "pdf")

## plot retention time for seed1 (all + mapped, pos and neg)
setwd("~/Documents/01_GWAS/01_Data/00_two_biological_replicates/")
rel_rep1_pos <- read.csv("rep1_positive_match.csv", sep="\t")
rel_rep1_neg <- read.csv("rep1_negative_match.csv", sep="\t")

rt_pos <- rel_rep1_pos[rel_rep1_pos[, "Name"] %in% rownames(heritability_pos), "RT"]
rt_pos_map <- rel_rep1_pos[rel_rep1_pos[, "Name"] %in% inds_map_pos, "RT"]
rt_neg <- rel_rep1_neg[rel_rep1_neg[, "Name"] %in% rownames(heritability_neg), "RT"]
rt_neg_map <- rel_rep1_neg[rel_rep1_neg[, "Name"] %in% inds_map_neg, "RT"]

df <- data.frame(
    values = c(rt_pos, rt_pos_map, rt_neg, rt_neg_map),
    ion_mode = c(rep("pos", length(c(rt_pos, rt_pos_map))),
        rep("neg", length(c(rt_neg, rt_neg_map)))),
    type = c(rep("all", length(rt_pos)), rep("mapped", length(rt_pos_map)),
        rep("all", length(rt_neg)), rep("mapped", length(rt_neg_map))))

g <- ggplot(df) + 
    facet_grid(ion_mode ~.) +
    geom_density(aes(x = values, color = type, fill = type), 
        alpha = .2, adjust = 1/20, position = "stack") +
    theme_bw() + xlab("retention time [min]") + 
    ylab("density (number of mass features)")

setwd("~/GitHub/GWAS_arabidopsis_seed/heritability")
ggsave(g, filename = "plot_density_rettime_mass_features.pdf", device = "pdf")

## truncate the table with mapping results (overlap) that it only contains those features
## meet the criteria (mz_dev, rt_dev, corr)
setwd("~/GitHub/GWAS_arabidopsis_seed/")
loci_pos <- read.table("gwas_complete_met_all_trueLociLOD_pos.txt", sep = "\t", header = TRUE)
loci_neg <- read.table("gwas_complete_met_all_trueLociLOD_neg.txt", sep = "\t", header = TRUE)

loci_pos_cut <- loci_pos[loci_pos[, "met_rep1"] %in% rownames(heritability_pos), ]
loci_neg_cut <- loci_neg[loci_neg[, "met_rep1"] %in% rownames(heritability_neg), ]

write.table(loci_pos_cut, "gwas_complete_met_all_trueLociLOD_pos_cut.txt", sep= "\t", quote = FALSE)
write.table(loci_neg_cut, "gwas_complete_met_all_trueLociLOD_neg_cut.txt", sep= "\t", quote = FALSE)
