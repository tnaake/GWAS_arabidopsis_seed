## load data sets
setwd("~/01_GWAS/00_MetNet")
peaklist_neg1 <- read.table("rep1_negative_match.csv", sep = "\t", dec = ".", stringsAsFactors = FALSE, header = TRUE)
peaklist_neg2 <- read.table("rep2_negative_match.csv", sep = "\t", dec = ".", stringsAsFactors = FALSE, header = TRUE)
peaklist_pos1 <- read.table("rep1_positive_match.csv", sep = "\t", dec = ".", stringsAsFactors = FALSE, header = TRUE)
peaklist_pos2 <- read.table("rep2_positive_match.csv", sep = "\t", dec = ".", stringsAsFactors = FALSE, header = TRUE)

## combine data sets (preparation)
pl_neg1 <- peaklist_neg1[, grep(colnames(peaklist_neg1), pattern = "Sample")]
colnames(pl_neg1) <- paste(colnames(pl_neg1), "_rep1", sep = "")
cols <- grep(colnames(peaklist_neg2), pattern = "Sample")
cols_l <- grep(colnames(peaklist_neg2), pattern = "leaf_Sample")
pl_neg2 <- peaklist_neg2[, cols[!cols %in% cols_l]]
colnames(pl_neg2) <- paste(colnames(pl_neg2), "_rep2", sep = "")

pl_pos1 <- peaklist_pos1[, grep(colnames(peaklist_pos1), pattern = "Sample")]
colnames(pl_pos1) <- paste(colnames(pl_pos1), "_rep1", sep = "")
cols <- grep(colnames(peaklist_pos2), pattern = "Sample")
cols_l <- grep(colnames(peaklist_pos2), pattern = "leaf_Sample")
pl_pos2 <- peaklist_pos2[, cols[!cols %in% cols_l]]
colnames(pl_pos2) <- paste(colnames(pl_pos2), "_rep2", sep = "")

## combine
pl_neg <- cbind(mz = peaklist_neg1[, "mz"], rt = peaklist_neg1[, "RT"], pl_neg1, pl_neg2)
pl_pos <- cbind(mz = peaklist_pos1[, "mz"], rt = peaklist_pos1[, "RT"], pl_pos1, pl_pos2)


## create data.frame with transformations
transformations <- rbind(
  c("transf_Hydrogenation/dehydrogenation", "H2", 2.0156500642, "?"), ## --> +
  c("transf_Acetylation (-H)", "C2H3O2", 59.0133043405, "?"), ## log P=-0.17 (acetic acid, Pubchem) --> +
  c("transf_Acetylation (-H2O)", "C2H2O",  42.0105646863, "?"), ## log P=-0.17 (acetic acid, Pubchem) --> +
  c("transf_benzoyl", "C6H4CO", 104.026213, "?"), ## --> +
  c("transf_hydroxy benzoyl", "C7H4O2", 120.021128, "?"), ## --> +
  c("transf_hydroxy benzyl", "C7H6O", 106.041863, "?"), ## --> +
  c("transf_galloyl", "C7H4O4", 152.010958, "?"), ## --> +
  c("transf_methylation", "CH2", 14.015649, "?"), ## --> + 
  c("transf_methoxylation", "CH2O", 30.010564, "?"),
  c("transf_CHO2", "CHO2", 44.9976542763, "?"),
  c("transf_Glyoxylate (-H2O)", "C2O2",  55.9898292442, "?"), ## log Kow=-0.07 (glyoxylic acid, Pubchem)
  c("transf_Biotinyl (-H)", "C10H15N2O3S", 243.0803380482, "?"), ## log Kow=0.39 (biotin, Pubchem) --> +
  c("transf_Biotinyl (-H2O)", "C10H14N2O2S", 226.0775983940, "?"), ## log Kow=0.39 (biotin, Pubchem) --> +
  c("transf_C2H2", "C2H2", 26.0156500642, "?"),
  c("transf_Sulfate (-H2O)", "SO3", 79.9568145563, "?"), 
  c("transf_Isoprene addition (-H)", "C5H7", 67.0547752247, "+"), ## log Kow=2.42 (isoprene, Pubchem)
  c("transf_Ketol group (-H2O)", "C2H2O", 42.0105646863, "?"),
  c("transf_Primary amine", "NH2", 16.0187240694, "?"),
  c("transf_Secondary amine", "NH", 15.0108990373, "?"),
  c("transf_Tertiary amine", "N", 14.0030740052, "?"),
  c("transf_Hydroxylation (-H)", "O", 15.9949146221, "-"),
  c("transf_Malonyl group (-H2O)", "C3H2O3", 86.0003939305, "?"), ## log Kow=-0.81 (malonic acid, Pubchem) --> -
  c("transf_Urea addition (-H)", "CH3N2O", 59.0245377288, "-"), ## log Kow=-2.11 (urea, Pubchem)
  c("transf_D-ribose (-H2O) (ribosylation)", "C5H8O4", 132.0422587452, "-"),
  c("transf_Rhamnose (-H2O)", "C6H10O4", 146.0579088094, "-"),
  c("transf_Disaccharide (-H2O) #1", "C12H20O10", 324.105649, "-"),
  c("transf_Disaccharide (-H2O) #2", "C12H20O11", 340.1005614851, "-"),  
  c("transf_Glucuronic acid (-H2O)", "C6H8O6", 176.0320879894, "-"), ## log P=-2.57 (glucoronic acid, Pubchem)
  c("transf_Monosaccharide (-H2O)", "C6H10O5", 162.0528234315, "-"), ## log Kow=-3.24 (glucose, Pubchem)
  c("transf_Deoxyhexose (-H20)", "C6H10O4", 146.0579088094, "-"),
  c("transf_Pentose (-H2O)", "C5H8O4", 132.042260, "-"),
  c("transf_Trisaccharide (-H2O)", "C18H30O15", 486.1584702945, "-"),
  c("transf_Glucose-O-phosphate (-H2O)", "C6H11O8P", 242.0191538399, "-"), 
  c("transf_coumaroyl (-H2O)", "C9H6O2", 146.0367794368, "+"), ## log P=1.79 (coumaric acid, Pubchem)
  c("transf_coumaroyl hexose", "C15H16O7", 308.089599, "?"), ## --> -
  c("transf_caffeoyl", "C9H6O3", 162.031693, "+"),
  c("transf_feruloyl (-H2O)", "C9H6O2OCH2", 176.0473441231, "+"), ## log P=1.51 (ferulic acid, Pubchem)
  c("transf_sinapoyl (-H2O)", "C9H6O2OCH2OCH2", 206.0579088094, "+"), 
  c("transf_protocatechuoyl", "C7H4O3", 136.016043, "?"), ## --> +
  c("transf_quinic acid (-H2O)", "C7H10O5", 174.052824, "?"), ## logP=-2.007 (quinic acid, chemspicer) --> -
  c("transf_shikimic acid (-H2O)", "C7H8O4", 156.042260, "?"), 
  c("transf_ellagic acid (-H2O)", "	C14H4O7", 283.995705, "?")) ## log Kow=-2.05 (ellagic acid, Pubchem) --> -
##c("putrescine to spermidine (+C3H7N)", "C3H7N", 57.0578492299, "?"))

## check within each pcgroup, mass differences with respect to the M+H
adducts_pos <- rbind(
  c("adduct_formic acid adduct", "CH2O2",  46.0054792, "?"), ## only positive
  c("adduct_ammonium", "NH4-H", 18.034374-1.007825, "?"),
  c("adduct_acetonitril H+CH3CN", "CH3CN", 41.026549, "?"),
  c("adduct_acetonitril H+CH3CN+1H2O", "CH3CNH2O", 59.037113, "?"),
  c("adduct_sodium formate adduct", "NaHCO2", 67.9874246, "?"),
  c("adduct_Na adduct (+Na+)", "Na-H", 22.989770-1.007825, "?"), ## only positive
  c("adduct_K adduct (+K+)", "K-H", 38.963708-1.007825, "?"), ## only positive
  c("adduct_isotopic+1_H", "isotopic peak 2H", 1.00628, "?"),
  c("adduct_isotopic+2_15N-14N", "isotopic peak 15N-14N", 0.997035, "?"),
  c("adduct_isotopic+2_18O-16O", "isotopic peak 18O-16O", 2.004244, "?"),
  c("adduct_isotopic+2_34S-32S", "isotopic peak 34S-32S", 1.995796, "?"),
  c("adduct_isotopic+1", "isotopic peak 13C1", 1.0033554, "?"),
  c("adduct_isotopic+2", "isotopic peak 13C2", 2.0067108, "?"),
  c("adduct_isotopic+3", "isotopic peak 13C3", 3.0100662, "?"),
  c("adduct_loss of deoxyhexose (e.g. Rhamnose)", "C6H10O4", -146.0579088094, "?"),
  c("adduct_loss of hexose", "C5H10O5", -162.0528234315, "?"),
  c("adduct_loss of pentose", "C5H8O4", -132.042260, "?"),
  c("adduct_loss of malonyl group", "C3H2O3", -86.0003939305, "?"),
  c("adduct_decarboxylation (loss from malonyl)", "CO2", -43.989830, "?"),
  c("adduct_loss of water", "H2O", -18.010565, "?"), 
  c("adduct_loss of COCH2", "COCH2", -42.0105646863, "?")
)

## check within each pcgroup, mass differences with respect to the M-H
adducts_neg <- rbind(
  c("adduct_formate HCOO-", "HCOO+H", 44.997654+1.007825, "?"),
  c("adduct_chlorine adduct", "Cl+H", 34.968853+1.007825, "?"),
  c("adduct_M-2H+Na", "Na-H", 22.989770-1.00782, "?"),
  c("adduct_sodium formate adduct", "NaHCO2", 67.9874246, "?"),
  c("adduct_isotopic+1_H", "isotopic peak 2H", 1.00628, "?"),
  c("adduct_isotopic+2_15N-14N", "isotopic peak 15N-14N", 0.997035, "?"),
  c("adduct_isotopic+2_18O-16O", "isotopic peak 18O-16O", 2.004244, "?"),
  c("adduct_isotopic+2_34S-32S", "isotopic peak 34S-32S", 1.995796, "?"),
  c("adduct_isotopic+1", "isotopic peak 13C1", 1.0033554, "?"),
  c("adduct_isotopic+2", "isotopic peak 13C2", 2.0067108, "?"),
  c("adduct_isotopic+3", "isotopic peak 13C3", 3.0100662, "?"), 
  c("adduct_loss of deoxyhexose (e.g. Rhamnose)", "C6H10O4", -146.0579088094, "?"), 
  c("adduct_loss of hexose", "C5H10O5", -162.0528234315, "?"),
  c("adduct_loss of pentose", "C5H8O4", -132.042260, "?"),
  c("adduct_loss of malonyl group", "C3H2O3", -86.0003939305, "?"),
  c("adduct_decarboxylation (loss from malonyl)", "CO2", -43.989830, "?"),
  c("adduct_loss of water", "H2O", -18.010565, "?"), 
  c("adduct_loss of COCH2", "COCH2", -42.0105646863, "?")
)

## positive mode
transformations_pos <- data.frame(
  group = c(adducts_pos[, 1], transformations[, 1]),
  formula = c(adducts_pos[, 2], transformations[, 2]), 
  mass = c(as.numeric(adducts_pos[, 3]), as.numeric(transformations[, 3])),
  rt = c(adducts_pos[, 4], transformations[, 4]))

## negative mode 
transformations_neg <- data.frame(
  group = c(adducts_neg[, 1], transformations[, 1]),
  formula = c(adducts_neg[, 2], transformations[, 2]), 
  mass = c(as.numeric(adducts_neg[, 3]), as.numeric(transformations[, 3])),
  rt = c(adducts_neg[, 4], transformations[, 4]))

## use function createStructuralAdjacency and remove false positives by 
## function rtCorrection
## pos
struct_adj_pos <- structural(x = pl_pos, 
    transformation = transformations_pos, ppm = 10, directed = TRUE)
struct_adj_pos <- rtCorrection(structural = struct_adj_pos, x = pl_pos, 
    transformation = transformations_pos)
## neg
struct_adj_neg <- structural(x = pl_neg, 
    transformation = transformations_neg, ppm = 10, directed = TRUE)
struct_adj_neg <- rtCorrection(structural = struct_adj_neg, x = pl_neg,
    transformation = transformations_neg)

rownames(struct_adj_neg[[1]]) <- colnames(struct_adj_neg[[1]]) <- peaklist_neg1[, "Name"]
rownames(struct_adj_neg[[2]]) <- colnames(struct_adj_neg[[2]]) <- peaklist_neg1[, "Name"]
rownames(struct_adj_pos[[1]]) <- colnames(struct_adj_pos[[1]]) <- peaklist_pos1[, "Name"]
rownames(struct_adj_pos[[2]]) <- colnames(struct_adj_pos[[2]]) <- peaklist_pos1[, "Name"]
## save
save(struct_adj_pos, file = "MetNet_seed_struct_adj_pos.RData")
save(struct_adj_neg, file = "MetNet_seed_struct_adj_neg.RData")

## use function statistical/threshold
inds_pos <- which(colnames(pl_pos) == "Sample1_Positive_rep1"):which(colnames(pl_pos) == "Sample338_Positive_rep2")
inds_neg <- which(colnames(pl_neg) == "Sample1_Negative_rep1"):which(colnames(pl_neg) == "Sample338_Negative_rep2")
pl_pos_cut <- pl_pos[, inds_pos]
pl_neg_cut <- pl_neg[, inds_neg]

models <- c("pearson", "pearson_partial", "spearman", "spearman_partial")

## apply the function statistical to create weighted adjacency matrices 
## per model
stat_adj_pos <- statistical(as.matrix(pl_pos_cut), model = models) 
stat_adj_neg <- statistical(as.matrix(pl_neg_cut), model = models)

rownames(stat_adj_pos[[1]]) <- colnames(stat_adj_pos[[1]]) <- peaklist_neg1[, "Name"]
rownames(stat_adj_pos[[2]]) <- colnames(stat_adj_pos[[2]]) <- peaklist_neg1[, "Name"]
rownames(stat_adj_pos[[3]]) <- colnames(stat_adj_pos[[3]]) <- peaklist_neg1[, "Name"]
rownames(stat_adj_pos[[4]]) <- colnames(stat_adj_pos[[4]]) <- peaklist_neg1[, "Name"]
rownames(stat_adj_neg[[1]]) <- colnames(stat_adj_neg[[1]]) <- peaklist_neg1[, "Name"]
rownames(stat_adj_neg[[2]]) <- colnames(stat_adj_neg[[2]]) <- peaklist_neg1[, "Name"]
rownames(stat_adj_neg[[3]]) <- colnames(stat_adj_neg[[3]]) <- peaklist_neg1[, "Name"]
rownames(stat_adj_neg[[4]]) <- colnames(stat_adj_neg[[4]]) <- peaklist_neg1[, "Name"]

save(stat_adj_pos, file = "MetNet_seed_stat_adj_pos.RData")
save(stat_adj_neg, file = "MetNet_seed_stat_adj_neg.RData")

## apply the function threshold to create unweighted adjacency matrices
## type = "threshold" 
## define thresholds
pdf("hist_pos_pearson.pdf")
hist(stat_adj_pos_NA[["pearson"]])
dev.off()

pdf("hist_pos_pearson_partial.pdf")
hist(stat_adj_pos_NA[["pearson_partial"]])
dev.off()

pdf("hist_pos_spearman.pdf")
hist(stat_adj_pos_NA[["spearman"]])
dev.off()

pdf("hist_pos_spearman_partial.pdf")
hist(stat_adj_pos_NA[["spearman_partial"]])
dev.off()

pdf("hist_neg_pearson.pdf")
hist(stat_adj_neg_NA[["pearson"]])
dev.off()

pdf("hist_neg_pearson_partial.pdf")
hist(stat_adj_neg_NA[["pearson_partial"]])
dev.off()

pdf("hist_neg_spearman.pdf")
hist(stat_adj_neg_NA[["spearman"]])
dev.off()

pdf("hist_neg_spearman_partial.pdf")
hist(stat_adj_neg_NA[["spearman_partial"]])
dev.off()

## check thresholds
table(stat_adj_neg[["pearson"]][upper.tri(stat_adj_neg[["pearson"]])] > 0.5)
table(stat_adj_pos[["pearson"]][upper.tri(stat_adj_pos[["pearson"]])] > 0.5)
table(stat_adj_neg[["pearson_partial"]][upper.tri(stat_adj_neg[["pearson_partial"]])] > 0.2)
table(stat_adj_pos[["pearson_partial"]][upper.tri(stat_adj_pos[["pearson_partial"]])] > 0.2)
table(stat_adj_neg[["spearman"]][upper.tri(stat_adj_neg[["spearman"]])] > 0.5)
table(stat_adj_pos[["spearman"]][upper.tri(stat_adj_pos[["spearman"]])] > 0.5)
table(stat_adj_neg[["spearman_partial"]][upper.tri(stat_adj_neg[["spearman_partial"]])] > 0.2)
table(stat_adj_pos[["spearman_partial"]][upper.tri(stat_adj_pos[["spearman_partial"]])] > 0.2)

args <- list("pearson" = 0.5, "pearson_partial" = 0.2, "spearman" = 0.5,
             "spearman_partial" = 0.2, threshold = 1)
stat_adj_pos_thr <- threshold(statistical = stat_adj_pos, type = "threshold",
                              args = args)
stat_adj_neg_thr <- threshold(statistical = stat_adj_neg, type = "threshold",
                              args = args)

## type = "top2"
args_top <- list(n = 250000)
stat_adj_pos_top2 <- threshold(statistical = stat_adj_pos, type = "top2",
                               args = args_top)
stat_adj_neg_top2 <- threshold(statistical = stat_adj_neg, type = "top2",
                               args = args_top)

save(stat_adj_pos_thr, stat_adj_pos_top2, file = "MetNet_seed_stat_adj_thr_pos.RData")
save(stat_adj_neg_thr, stat_adj_neg_top2, file = "MetNet_seed_stat_adj_thr_neg.RData")

## use function combine to combine the structural and statistical information
cons_adj_pos <- combine(structural = struct_adj_pos, 
                        statistical = stat_adj_pos_thr)
cons_adj_neg <- combine(structural = struct_adj_neg, 
                        statistical = stat_adj_neg_thr)

save(cons_adj_pos, file = "MetNet_strawberry_cons_adj_pos.RData")
save(cons_adj_neg, file = "MetNet_strawberry_cons_adj_neg.RData")

## rename edges, Cluster name from rep1
rownames(cons_adj_pos[[1]]) <- colnames(cons_adj_pos[[1]]) <- paste(peaklist_pos1[, "Name"], round(peaklist_pos1[, "mz"], 4), round(peaklist_pos1[, "RT"], 2), sep = "_") 
rownames(cons_adj_pos[[2]]) <- colnames(cons_adj_pos[[2]]) <- paste(peaklist_pos1[, "Name"], round(peaklist_pos1[, "mz"], 4), round(peaklist_pos1[, "RT"], 2), sep = "_")
rownames(cons_adj_neg[[1]]) <- colnames(cons_adj_neg[[1]]) <- paste(peaklist_neg1[, "Name"], round(peaklist_neg1[, "mz"], 4), round(peaklist_neg1[, "RT"], 2), sep = "_")
rownames(cons_adj_neg[[2]]) <- colnames(cons_adj_neg[[2]]) <- paste(peaklist_neg1[, "Name"], round(peaklist_neg1[, "mz"], 4), round(peaklist_neg1[, "RT"], 2), sep = "_")

## remove adducts
cons_adj_pos_rem <- cons_adj_pos
cons_adj_neg_rem <- cons_adj_neg
rown <- rownames(cons_adj_pos[[1]])
rt <- as.numeric(unlist(lapply(strsplit(rown, split = "_"), "[", 4)))
for (i in 1:nrow(cons_adj_pos_rem[[1]])) {
  cons_adj_pos_rem[[1]][abs(rt[i] - rt) > 0.1 & grepl(cons_adj_pos_rem[[1]][, i], pattern = "adduct_"), ] <- 0
  cons_adj_pos_rem[[2]][abs(rt[i] - rt) > 0.1 & grepl(cons_adj_pos_rem[[2]][, i], pattern = "adduct_"), ] <- 0
}
rown <- rownames(cons_adj_neg[[1]])
rt <- as.numeric(unlist(lapply(strsplit(rown, split = "_"), "[", 4)))
for (i in 1:nrow(cons_adj_neg_rem[[1]])) {
    cons_adj_neg_rem[[1]][abs(rt - rt[i]) > 0.1 & grepl(cons_adj_neg_rem[[2]][, i], pattern = "adduct_"), ] <- 0
    cons_adj_neg_rem[[2]][abs(rt - rt[i]) > 0.1 & grepl(cons_adj_neg_rem[[2]][, i], pattern = "adduct_"), ] <- ""
}

## build network with mapped features 
## (feature that was mapped at least to one locus in one of the replicates)
options(stringsAsFactors = FALSE)
map_pos <- read.table("gwas_complete_met_all_trueLociLOD_pos.txt", 
    sep = "\t", header = TRUE)
map_pos <- map_pos[map_pos[, "locus_tag_seed1"] != "" | map_pos[, "locus_tag_seed2"] != "", ]
map_pos_met <- unique(map_pos[, "met_rep1"])
map_neg <- read.table("gwas_complete_met_all_trueLociLOD_neg.txt", 
    sep = "\t", header = TRUE)
map_neg <- map_neg[map_neg[, "locus_tag_seed1"] != "" | map_neg[, "locus_tag_seed2"] != "", ]
map_neg_met <- unique(map_neg[, "met_rep1"])

## truncate that it only contains mapped features
cons_adj_pos_map <- cons_adj_pos_rem
cut_rown <- unlist(lapply(lapply(strsplit(rownames(cons_adj_pos_rem[[1]]), 
    split = "_"), "[", 1:2), paste, collapse = "_"))
inds <- match(cut_rown, map_pos_met)
inds <- !is.na(inds)
cons_adj_pos_map[[1]] <- cons_adj_pos[[1]][inds, inds]
cons_adj_pos_map[[2]] <- cons_adj_pos[[2]][inds, inds]

cons_adj_neg_map <- cons_adj_neg_rem
cut_rown <- unlist(lapply(lapply(strsplit(rownames(cons_adj_neg_rem[[1]]), 
    split = "_"), "[", 1:2), paste, collapse = "_"))
inds <- match(cut_rown, map_neg_met)
inds <- !is.na(inds)
cons_adj_neg_map[[1]] <- cons_adj_neg_rem[[1]][inds, inds]
cons_adj_neg_map[[2]] <- cons_adj_neg_rem[[2]][inds, inds]

g_pos <- graph_from_adjacency_matrix(cons_adj_pos_map[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_map[[1]], mode = "directed")

## export to XML
write_graph(g_pos, file = "cons_adj_pos_map_graphml.xml", format = "graphml")
write_graph(g_neg, file = "cons_adj_neg_map_graphml.xml", format = "graphml")
