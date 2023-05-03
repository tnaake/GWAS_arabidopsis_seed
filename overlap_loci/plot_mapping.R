library(ggplot2)
library(dplyr)
setwd("~/GitHub/GWAS_arabidopsis_seed/overlap_loci")

################################################################################
## negative

## load data sets
options(stringsAsFactors = FALSE)
map_neg <- read.table(
    "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_neg_cut.txt", 
    header = TRUE, sep = "\t")

## load the relation to rep1, rep2, rep_leaf: there is a link between rep2 and 
## rep_feng (column mapping_rep2_rep_feng)
rel_rep1 <- read.csv("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/rep1_negative_match.csv", sep = "\t")
rel_rep2 <- read.csv("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/rep2_negative_match.csv", sep = "\t")
rel_rep_feng <- read.csv("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/rep_feng_negative_match.csv", sep = "\t")

################################################################################
## Correlation plots of LOD values
## seed1/seed2
tmp <- map_neg |>
    filter(!(locus_tag_seed1 == "" | is.na(locus_tag_seed1))) |>
    filter(!(locus_tag_seed2 == "" | is.na(locus_tag_seed2)))
x <- tmp[, "bestSNP_lod_seed1"]
y <- tmp[, "bestSNP_lod_seed2"]
#ks.test(x, rnorm(length(x), mean(x), sd(x)))
#ks.test(y, rnorm(length(y), mean(y), sd(y)))
#shapiro.test(x)
#shapiro.test(y)
cor_test <- cor.test(x, y, method = "spearman")
p_s1s2_neg <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits = c(5.0, 40)) + 
    scale_y_continuous(limits = c(5.0, 40)) + 
    geom_rug(col = rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("highest LOD, seed (repl. 1)") + ylab("highest LOD, seed (repl. 2)") +
    ggtitle(bquote(rho == .(paste(round(cor_test$estimate[[1]], 3))))) +
    coord_fixed() + theme_bw() +
    theme(legend.position = "none", axis.title = element_text(size = 8), 
          plot.title = element_text(size = 10))
ggsave(p_s1s2_neg, file = "plot_neg_seed1_seed2_scatter.pdf", device = "pdf", 
    width = 16, height = 16.8, units = "cm", useDingbats = FALSE)

## seed1/leaf2
tmp <- map_neg |>
    filter(!(locus_tag_seed1 == "" | is.na(locus_tag_seed1))) |>
    filter(!(locus_tag_leaf2 == "" | is.na(locus_tag_leaf2)))
x <- tmp[, "bestSNP_lod_seed1"]
y <- tmp[, "bestSNP_lod_leaf2"]
#ks.test(x, rnorm(length(x), mean(x), sd(x)))
#ks.test(y, rnorm(length(y), mean(y), sd(y)))
#shapiro.test(x)
#shapiro.test(y)
cor_test <- cor.test(x, y, method = "spearman")
p_s1l2_neg <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits = c(5.0, 40)) + 
    scale_y_continuous(limits = c(5.0, 40)) + 
    geom_rug(col = rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("highest LOD, seed (repl. 1)") + ylab("highest LOD, leaf") +
    ggtitle(bquote(rho == .(paste(round(cor_test$estimate[[1]], 3))))) +
    coord_fixed() + theme_bw() +
    theme(legend.position = "none", axis.title = element_text(size = 8), 
        plot.title = element_text(size = 10))
ggsave(p_s1l2_neg, file = "plot_neg_seed1_leaf2_scatter.pdf", device = "pdf", 
    width = 16, height = 16.8, units = "cm", useDingbats = FALSE)

## seed1/leaf_zhu
tmp <- map_neg |>
    filter(!(locus_tag_seed1 == "" | is.na(locus_tag_seed1))) |>
    filter(!(locus_tag_leaf_feng == "" | is.na(locus_tag_leaf_feng)))
x <- tmp[, "bestSNP_lod_seed1"]
y <- tmp[, "bestSNP_lod_leaf_feng"]
#ks.test(x, rnorm(length(x), mean(x), sd(x)))
#ks.test(y, rnorm(length(y), mean(y), sd(y)))
#shapiro.test(x)
#shapiro.test(y)
cor_test <- cor.test(x, y, method = "spearman")
p_s1lz_neg <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits = c(5.0, 40)) + 
    scale_y_continuous(limits = c(5.0, 40)) + 
    geom_rug(col = rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("highest LOD, seed (repl. 1)") + ylab("highest LOD, leaf \n(Zhu et al., 2022)") +
    ggtitle(bquote(rho == .(paste(round(cor_test$estimate[[1]], 3))))) +
    coord_fixed() + theme_bw() + 
    theme(legend.position = "none", axis.title = element_text(size = 8), 
        plot.title = element_text(size = 10))
ggsave(p_s1lz_neg, file = "plot_neg_seed1_leaf_zhu_scatter.pdf", device = "pdf", 
    width = 16, height = 16.8, units = "cm", useDingbats = FALSE)

## seed2/leaf2
tmp <- map_neg |>
    filter(!(locus_tag_seed2 == "" | is.na(locus_tag_seed2))) |>
    filter(!(locus_tag_leaf2 == "" | is.na(locus_tag_leaf2)))
x <- tmp[, "bestSNP_lod_seed2"]
y <- tmp[, "bestSNP_lod_leaf2"]
#ks.test(x, rnorm(length(x), mean(x), sd(x)))
#ks.test(y, rnorm(length(y), mean(y), sd(y)))
#shapiro.test(x)
#shapiro.test(y)
cor_test <- cor.test(x, y, method = "spearman")
p_s2l2_neg <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits = c(5.0, 40)) + 
    scale_y_continuous(limits = c(5.0, 40)) + 
    geom_rug(col = rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("highest LOD, seed (repl. 2)") + ylab("highest LOD, leaf") +
    ggtitle(bquote(rho == .(paste(round(cor_test$estimate[[1]], 3))))) +
    coord_fixed() + theme_bw() +
    theme(legend.position = "none", axis.title = element_text(size = 8), 
          plot.title = element_text(size = 10))
ggsave(p_s2l2_neg, file = "plot_neg_seed2_leaf2_scatter.pdf", device = "pdf", 
    width = 16, height = 16.8, units = "cm", useDingbats = FALSE)

## seed2/leaf_zhu
tmp <- map_neg  |>
    filter(!(locus_tag_seed2 == "" | is.na(locus_tag_seed2))) |>
    filter(!(locus_tag_leaf_feng == "" | is.na(locus_tag_leaf_feng)))
x <- tmp[, "bestSNP_lod_seed2"]
y <- tmp[, "bestSNP_lod_leaf_feng"]
#ks.test(x, rnorm(length(x), mean(x), sd(x)))
#ks.test(y, rnorm(length(y), mean(y), sd(y)))
#shapiro.test(x)
#shapiro.test(y)
cor_test <- cor.test(x, y, method = "spearman")
p_s2lz_neg <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits = c(5.0, 40)) + 
    scale_y_continuous(limits = c(5.0, 40)) + 
    geom_rug(col = rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("highest LOD, seed (repl. 2)") + ylab("highest LOD, leaf \n(Zhu et al., 2022)") +
    ggtitle(bquote(rho == .(paste(round(cor_test$estimate[[1]], 3))))) +
    coord_fixed() + theme_bw() + 
    theme(legend.position = "none", axis.title = element_text(size = 8), 
        plot.title = element_text(size = 10))
ggsave(p_s2lz_neg, file = "plot_neg_seed2_leaf_zhu_scatter.pdf", device = "pdf", 
     width = 16, height = 16.8, units = "cm", useDingbats = FALSE)

## seed2/leaf2
tmp <- map_neg  |>
    filter(!(locus_tag_seed2 == "" | is.na(locus_tag_seed2))) |>
    filter(!(locus_tag_leaf2 == "" | is.na(locus_tag_leaf2)))
x <- tmp[, "bestSNP_lod_seed2"]
y <- tmp[, "bestSNP_lod_leaf2"]
#ks.test(x, rnorm(length(x), mean(x), sd(x)))
#ks.test(y, rnorm(length(y), mean(y), sd(y)))
#shapiro.test(x)
#shapiro.test(y)
cor_test <- cor.test(x, y, method = "spearman")
p_s2l2_neg <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits = c(5.0, 40)) + 
    scale_y_continuous(limits = c(5.0, 40)) + 
    geom_rug(col = rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("highest LOD, seed (repl. 2)") + ylab("highest LOD, leaf") +
    ggtitle(bquote(rho == .(paste(round(cor_test$estimate[[1]], 3))))) +
    coord_fixed() + theme_bw() + 
    theme(legend.position = "none", axis.title = element_text(size = 8), 
        plot.title = element_text(size = 10))
ggsave(p_s2l2_neg, file = "plot_neg_seed2_leaf2_scatter.pdf", device = "pdf", 
    width = 16, height = 16.8, units = "cm", useDingbats = FALSE)

## leaf2/leaf_zhu
tmp <- map_neg  |>
    filter(!(locus_tag_leaf2 == "" | is.na(locus_tag_leaf2))) |>
    filter(!(locus_tag_leaf_feng == "" | is.na(locus_tag_leaf_feng)))
x <- tmp[, "bestSNP_lod_leaf2"]
y <- tmp[, "bestSNP_lod_leaf_feng"]
#ks.test(x, rnorm(length(x), mean(x), sd(x)))
#ks.test(y, rnorm(length(y), mean(y), sd(y)))
#shapiro.test(x)
#shapiro.test(y)
cor_test <- cor.test(x, y, method = "spearman")
p_l2lz_neg <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits = c(5.0, 40)) + 
    scale_y_continuous(limits = c(5.0, 40)) + 
    geom_rug(col = rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("highest LOD, leaf") + ylab("highest LOD, leaf \n(Zhu et al., 2022)") +
    ggtitle(bquote(rho == .(paste(round(cor_test$estimate[[1]], 3))))) +
    coord_fixed() + theme_bw() +
    theme(legend.position = "none", axis.title = element_text(size = 8), 
        plot.title = element_text(size = 10))
ggsave(p_l2lz_neg, file = "plot_neg_leaf2_leaf_zhu_scatter.pdf", device = "pdf", 
    width = 16, height = 16.8, units = "cm", useDingbats = FALSE)



################################################################################
## UpSetR
library(UpSetR)

## negative
## the map_neg will now be truncated such that they only contain the mass
## features of the core set
mapping_rep1_rep2 <- strsplit(rel_rep2[, "mapping_rep1_rep2"], split = "/")
mapping_df <- data.frame(
    rep1 = unlist(lapply(mapping_rep1_rep2, "[", 1)),
    rep2 = unlist(lapply(mapping_rep1_rep2, "[", 2))
)
mapping_rep2_rep_feng <- strsplit(rel_rep_feng[, "mapping_rep2_rep_feng"], 
                                  split = "/")
mapping_df_tmp <- data.frame(
    rep2 = unlist(lapply(mapping_rep2_rep_feng,"[", 1)),
    rep_feng = unlist(lapply(mapping_rep2_rep_feng, "[", 2))
)
mapping_df_all <- dplyr::inner_join(mapping_df, mapping_df_tmp, 
                                    by = "rep2", multiple = "all")
mapping_df_all <- mapping_df_all[!duplicated(
    paste(mapping_df_all$rep1, mapping_df_all$rep2, mapping_df_all$rep_feng)), ]

## create the UpSet plots

## first create binary matrices, fill the entries if 1 to the number of 
## row to mimic unique names

## original file
# cols <- c("locus_tag_seed1", "locus_tag_seed2", "locus_tag_leaf2", 
#     "locus_tag_leaf_feng", "met_rep1", "met_rep2")
# binary_mat <- map_neg |>
#     dplyr::filter(met_rep1 %in% mapping_df$rep1) |>
#     dplyr::filter(met_rep2 %in% mapping_df$rep2) |>
#     
#     dplyr::select(all_of(cols))
# binary_mat[, 1:4] <- ifelse(is.na(binary_mat[, 1:4]) | binary_mat[, 1:4] == "", 0, 1) |>
#     as.data.frame()

## for normalized
cols <- c("locus_tag_seed1", "locus_tag_seed2", "locus_tag_leaf2",
          "locus_tag_leaf_feng", "met_rep1", "met_rep2", "met_leaf_feng")
binary_mat <- map_neg |>
    dplyr::filter(
        met_rep1 %in% mapping_df_all$rep1 | is.na(met_rep1)) |>
    dplyr::filter(
        met_rep2 %in% mapping_df_all$rep2 | is.na(met_rep2)) |>
    dplyr::filter(
        met_leaf_feng %in% mapping_df_all$rep_feng | is.na(met_leaf_feng)) |>
    dplyr::select(all_of(cols))
binary_mat[, 1:4] <- ifelse(
    is.na(binary_mat[, 1:4]) | binary_mat[, 1:4] == "", 0, 1) |>
    as.data.frame()
# 
# ## for batch
# cols <- c("locus_tag_seed1", "locus_tag_seed2", "locus_tag_leaf2", 
#           "locus_tag_leaf_feng", "met_rep1", "met_rep2")
# binary_mat_batch <- trueLociLOD_batch |>
#     dplyr::filter(
#         met_rep1 %in% mapping_df_all$rep1 | is.na(met_rep1)) |>
#     dplyr::filter(
#         met_rep2 %in% mapping_df_all$rep2 | is.na(met_rep2)) |>
#     dplyr::filter(
#         met_leaf_feng %in% mapping_df_all$rep_feng | is.na(met_leaf_feng)) |>
#     dplyr::select(all_of(cols))
# binary_mat_batch[, 1:4] <- ifelse(
#     is.na(binary_mat_batch[, 1:4]) | binary_mat_batch[, 1:4] == "", 0, 1) |>
#     as.data.frame()

## do the actual plotting
binary_mat_upset <- binary_mat
colnames(binary_mat_upset)[1:4] <- c("seed \n (repl. 1)", "seed \n (repl. 2)", "leaf", "leaf \n (Zhu et al., 2022)")
m <- ComplexHeatmap::make_comb_mat(binary_mat_upset[, 1:4], mode = "distinct")

p_upset_neg <- ComplexHeatmap::UpSet(m = m, 
    comb_order = order(ComplexHeatmap::comb_size(m), decreasing = TRUE),
    top_annotation = upset_top_annotation(m, add_numbers = TRUE, numbers_rot = 90, 
        height = unit(13, "cm"), annotation_name_rot = 90),
    right_annotation = upset_right_annotation(m, add_numbers = TRUE),
    column_names_max_height = unit(1, "cm"))


## arrange in one plot
p_1 <- ggarrange(p_s1s2_neg, p_l2lz_neg, p_s2l2_neg, 
    p_s1lz_neg, p_s1l2_neg, p_s2lz_neg,
    ncol = 3, nrow = 2, labels = "AUTO")
p_2 <- ggarrange(grid::grid.grabExpr(ComplexHeatmap::draw(p_upset_neg)), ncol = 1, nrow = 1, labels = "H")

p <- ggarrange(p_1, p_2, ncol = 1, nrow = 2, 
               labels = NULL, heights = c(1.4, 2))
ggsave(p, filename = "figure_mapping_intersection_neg.pdf", width = 210, height = 297, units = "mm")

################################################################################
plot_loci_distribution <- function(binary_mat = binary_mat,
    s1 = FALSE, s2 = FALSE, l2 = FALSE, lz = FALSE, file = "plot.pdf") {
    
    s1 <- ifelse(s1, 1, 0)
    s2 <- ifelse(s2, 1, 0)
    l2 <- ifelse(l2, 1, 0)
    
    if (!is.null(lz)) {
        lz <- ifelse(lz, 1, 0)
        binary_mat$met <- paste(binary_mat$met_rep1, binary_mat$met_rep2, 
                                binary_mat$met_leaf_feng)
        df <- data.frame(values = table(binary_mat[
            which(binary_mat[, "locus_tag_seed1"] == s1 & 
                      binary_mat[, "locus_tag_seed2"] == s2 & 
                      binary_mat[, "locus_tag_leaf2"] == l2 &
                      binary_mat[, "locus_tag_leaf_feng"] == lz), "met"]))
    } else {
        binary_mat$met <- paste(binary_mat$met_rep1, binary_mat$met_rep2)
        df <- data.frame(values = table(
            binary_mat[which(binary_mat[, "locus_tag_seed1"] == s1 & 
                binary_mat[, "locus_tag_seed2"] == s2 & 
                binary_mat[, "locus_tag_leaf2"] == l2), "met"]))
    }
    p <- ggplot(df, aes(x = values.Freq)) + 
        geom_bar(aes(y = (after_stat(count)) / sum(after_stat(count)))) + 
        ylim(c(0, 1)) +xlim(c(0, 75)) + 
        xlab("") + ylab("") +
        theme_bw()
    ggsave(p, file = file, device = "pdf", useDingbats = FALSE)
    p
}

## only s1
p_s1_neg <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = FALSE, l2 = FALSE, lz = FALSE, 
    file = "plot_distribution_neg_number_mass_features_s1.pdf")

## only s2
p_s2_neg <- plot_loci_distribution(binary_mat, 
    s1 = FALSE, s2 = TRUE, l2 = FALSE, lz = FALSE, 
    file = "plot_distribution_neg_number_mass_features_s2.pdf")

## only l2
p_l2_neg <- plot_loci_distribution(binary_mat, 
    s1 = FALSE, s2 = FALSE, l2 = TRUE, lz = FALSE, 
    file = "plot_distribution_neg_number_mass_features_l2.pdf")

## only leaf_feng
p_lz_neg <- plot_loci_distribution(binary_mat, 
    s1 = FALSE, s2 = FALSE, l2 = FALSE, lz = TRUE, 
    file = "plot_distribution_neg_number_mass_features_s1.pdf")

## seed1&seed2
p_s1s2_neg <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = TRUE, l2 = FALSE, lz = FALSE, 
    file = "plot_distribution_neg_number_mass_features_s1s2.pdf")

## seed1&leaf2
p_s1l2_neg <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = FALSE, l2 = TRUE, lz = FALSE, 
    file = "plot_distribution_neg_number_mass_features_s1l2.pdf")

## seed1&leaf_feng
p_s1lz_neg <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = FALSE, l2 = FALSE, lz = TRUE, 
    file = "plot_distribution_neg_number_mass_features_s1lz.pdf")

## seed2&leaf2
p_s2l2_neg <- plot_loci_distribution(binary_mat, 
    s1 = FALSE, s2 = TRUE, l2 = TRUE, lz = FALSE, 
    file = "plot_distribution_neg_number_mass_features_s2l2.pdf")

## seed2&leaf_feng
p_s2lz_neg <- plot_loci_distribution(binary_mat, 
    s1 = FALSE, s2 = TRUE, l2 = FALSE, lz = TRUE, 
    file = "plot_distribution_neg_number_mass_features_s2lz.pdf")

## seed1&seed2&leaf2
p_s1s2l2_neg <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = TRUE, l2 = TRUE, lz = FALSE, 
    file = "plot_distribution_neg_number_mass_features_s1s2l2.pdf")

## seed1&seed2&leaf_feng
p_s1s2lz_neg <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = TRUE, l2 = FALSE, lz = TRUE, 
    file = "plot_distribution_neg_number_mass_features_s1s2lz.pdf")

## seed1&leaf2&leaf_feng
p_s1l2lz_neg <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = FALSE, l2 = TRUE, lz = TRUE, 
    file = "plot_distribution_neg_number_mass_features_s1l2lz.pdf")

## seed2&leaf2&leaf_feng
p_s2l2lz_neg <- plot_loci_distribution(binary_mat, 
    s1 = FALSE, s2 = TRUE, l2 = TRUE, lz = TRUE, 
    file = "plot_distribution_neg_number_mass_features_s2l2lz.pdf")

## seed1&seed2&leaf2&leaf_feng
p_s1s2l2lz_neg <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = TRUE, l2 = TRUE, lz = TRUE, 
    file = "plot_distribution_neg_number_mass_features_s1s2l2lz.pdf")


################################################################################
## positive mode 
## load data sets
options(stringsAsFactors = FALSE)
map_pos <- read.table(
    "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_pos_cut.txt", 
    header = TRUE, sep = "\t")

## load the relation to rep1, rep2: there is a link between rep1 and 
## rep2 (column mapping_rep1_rep2)
rel_rep1 <- read.csv("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/rep1_positive_match.csv", sep = "\t")
rel_rep2 <- read.csv("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/rep2_positive_match.csv", sep = "\t")


################################################################################
## Correlation plots of LOD values

## seed1/seed2
tmp <- map_pos  |>
    filter(!(locus_tag_seed1 == "" | is.na(locus_tag_seed1))) |>
    filter(!(locus_tag_seed2 == "" | is.na(locus_tag_seed2)))
x <- tmp[, "bestSNP_lod_seed1"]
y <- tmp[, "bestSNP_lod_seed2"]
#ks.test(x, rnorm(length(x), mean(x), sd(x)))
#ks.test(y, rnorm(length(y), mean(y), sd(y)))
#shapiro.test(x)
#shapiro.test(y)

cor_test <- cor.test(x, y, method = "spearman")
p_s1s2_pos <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits = c(5.0, 40)) + 
    scale_y_continuous(limits = c(5.0, 40)) + 
    geom_rug(col = rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("highest LOD, seed (repl. 1)") + ylab("highest LOD, seed (repl. 2)") +
    ggtitle(bquote(rho == .(paste(round(cor_test$estimate[[1]], 3))))) +
    coord_fixed() + theme_bw() + 
    theme(legend.position = "none", axis.title = element_text(size = 8), 
        plot.title = element_text(size = 10))
ggsave(p_s1s2_pos, file = "plot_pos_seed1_seed2_scatter.pdf", device = "pdf", 
        width = 16, height = 16.8, units = "cm", useDingbats = FALSE)

## seed1/leaf2
tmp <- map_pos |>
    filter(!(locus_tag_seed1 == "" | is.na(locus_tag_seed1))) |>
    filter(!(locus_tag_leaf2 == "" | is.na(locus_tag_leaf2)))
x <- tmp[, "bestSNP_lod_seed1"]
y <- tmp[, "bestSNP_lod_leaf2"]
#ks.test(x, rnorm(length(x), mean(x), sd(x)))
#ks.test(y, rnorm(length(y), mean(y), sd(y)))
#shapiro.test(x)
#shapiro.test(y)

cor_test <- cor.test(x, y, method = "spearman")
p_s1l2_pos <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits = c(5.0, 40)) + 
    scale_y_continuous(limits = c(5.0, 40)) + 
    geom_rug(col = rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("highest LOD, seed (repl. 1)") + ylab("highest LOD, leaf") +
    ggtitle(bquote(rho == .(paste(round(cor_test$estimate[[1]], 3))))) +
    coord_fixed() + theme_bw() + 
    theme(legend.position = "none", axis.title = element_text(size = 8), 
          plot.title = element_text(size = 10))
ggsave(p_s1l2_pos, file = "plot_pos_seed1_leaf2_scatter.pdf", device = "pdf", 
       width = 16, height = 16.8, units = "cm", useDingbats = FALSE)

## seed2/leaf2
tmp <- map_pos  |>
    filter(!(locus_tag_seed2 == "" | is.na(locus_tag_seed2))) |>
    filter(!(locus_tag_leaf2 == "" | is.na(locus_tag_leaf2)))
x <- tmp[, "bestSNP_lod_seed2"]
y <- tmp[, "bestSNP_lod_leaf2"]
#ks.test(x, rnorm(length(x), mean(x), sd(x)))
#ks.test(y, rnorm(length(y), mean(y), sd(y)))
#shapiro.test(x)
#shapiro.test(y)

cor_test <- cor.test(x, y, method = "spearman")
p_s2l2_pos <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits = c(5.0, 40)) + 
    scale_y_continuous(limits = c(5.0, 40)) + 
    geom_rug(col = rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("highest LOD, seed (repl. 2)") + ylab("highest LOD, leaf") +
    ggtitle(bquote(rho == .(paste(round(cor_test$estimate[[1]], 3))))) +
    coord_fixed() + theme_bw() + 
    theme(legend.position = "none", axis.title = element_text(size = 8), 
          plot.title = element_text(size = 10))
ggsave(p_s2l2_pos, file = "plot_pos_seed2_leaf2_scatter.pdf", device = "pdf", 
       width = 16, height = 16.8, units = "cm", useDingbats = FALSE)

################################################################################

## UpSetR
library(UpSetR)

## positive
## the trueLociLOD will now be truncated such that they only contain the mass
## features of the core set
mapping_rep1_rep2 <- strsplit(rel_rep2[, "mapping_rep1_rep2"], split = "/")
mapping_df_all <- data.frame(
    rep1 = unlist(lapply(mapping_rep1_rep2, "[", 1)),
    rep2 = unlist(lapply(mapping_rep1_rep2, "[", 2))
)

## create the UpSet plots

## first create binary matrices, fill the entries if 1 to the number of 
## row to mimic unique names

## create binary_mat
cols <- c("locus_tag_seed1", "locus_tag_seed2", "locus_tag_leaf2", 
    "met_rep1", "met_rep2")
binary_mat <- map_pos |>
    dplyr::filter(met_rep1 %in% mapping_df$rep1) |>
    dplyr::filter(met_rep2 %in% mapping_df$rep2) |>
    dplyr::select(all_of(cols))
binary_mat[, 1:3] <- ifelse(is.na(binary_mat[, 1:3]) | binary_mat[, 1:3] == "", 0, 1) |>
    as.data.frame()

## do the actual plotting
binary_mat_upset <- binary_mat
colnames(binary_mat_upset)[1:3] <- c("seed \n (repl. 1)", "seed \n (repl. 2)", "leaf")
m <- ComplexHeatmap::make_comb_mat(binary_mat_upset[, 1:3], mode = "distinct")

p_upset_pos <- ComplexHeatmap::UpSet(m = m, 
    comb_order = order(ComplexHeatmap::comb_size(m), decreasing = TRUE),
    top_annotation = upset_top_annotation(m, add_numbers = TRUE, numbers_rot = 90, 
        height = unit(13, "cm"), annotation_name_rot = 90),
    right_annotation = upset_right_annotation(m, add_numbers = TRUE),
    column_names_max_height = unit(1, "cm"))

## arrange in one plot
p_1 <- ggarrange(p_s1s2_pos, p_s1l2_pos, p_s2l2_pos,
    ncol = 3, nrow = 1, labels = "AUTO")
p_2 <- ggarrange(grid::grid.grabExpr(ComplexHeatmap::draw(p_upset_pos)), 
    ncol = 1, nrow = 1, labels = "D")

p <- ggarrange(p_1, p_2, ncol = 1, nrow = 2, labels = NULL, heights = c(1.0, 2))
ggsave(p, filename = "figure_mapping_intersection_pos.pdf", width = 210, height = 267, units = "mm")


################################################################################
## only s1
p_s1_pos <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = FALSE, l2 = FALSE, lz = NULL, 
    file = "plot_distribution_pos_number_mass_features_s1.pdf")

## only s2
p_s2_pos <- plot_loci_distribution(binary_mat, 
    s1 = FALSE, s2 = TRUE, l2 = FALSE, lz = NULL, 
    file = "plot_distribution_pos_number_mass_features_s2.pdf")

## only l2
p_l2_pos <- plot_loci_distribution(binary_mat, 
    s1 = FALSE, s2 = FALSE, l2 = TRUE, lz = NULL, 
    file = "plot_distribution_pos_number_mass_features_l2.pdf")

## seed1&seed2
p_s1s2_pos <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = TRUE, l2 = FALSE, lz = NULL, 
    file = "plot_distribution_pos_number_mass_features_s1s2.pdf")

## seed1&leaf2
p_s1l2_pos <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = FALSE, l2 = TRUE, lz = NULL, 
    file = "plot_distribution_pos_number_mass_features_s1l2.pdf")

## seed2&leaf2
p_s2l2_pos <- plot_loci_distribution(binary_mat, 
    s1 = FALSE, s2 = TRUE, l2 = TRUE, lz = NULL, 
    file = "plot_distribution_pos_number_mass_features_s2l2.pdf")

## seed1&seed2&leaf2
p_s1s2l2_pos <- plot_loci_distribution(binary_mat, 
    s1 = TRUE, s2 = TRUE, l2 = TRUE, lz = NULL, 
    file = "plot_distribution_pos_number_mass_features_s1s2l2.pdf")


## create a plot with negative and positive mode

p <- ggarrange(p_s1_neg, p_s2_neg, p_l2_neg, p_lz_neg, p_s1s2_neg, p_s1l2_neg,
          p_s1lz_neg, p_s2l2_neg, p_s2lz_neg, p_s1s2l2_neg, p_s1s2lz_neg,
          p_s1l2lz_neg, p_s2l2lz_neg, p_s1s2l2lz_neg, 
          
          p_s1_pos, p_s2_pos, p_l2_pos, p_s1s2_pos, p_s1l2_pos, p_s2l2_pos,
          p_s1s2l2_pos,
          nrow = 7, ncol = 3, label = "AUTO")



