library(ggplot2)

## load data sets
options(stringsAsFactors = FALSE)
map_neg <- read.table(
    "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_neg_cut.txt", 
    header = TRUE, sep = "\t")
map_pos <- read.table(
    "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_pos_cut.txt", 
    header = TRUE, sep = "\t")

setwd("~/GitHub/GWAS_arabidopsis_seed/overlap_loci")

## negative
## seed1/seed2
map_tmp <- map_neg[!is.na(map_neg[, "locus_tag_seed1"]), ]
s1s2_neg <- map_tmp[map_tmp[, "locus_tag_seed1"] != "" & 
    map_tmp[, "locus_tag_seed2"] != "", ]
s1s2_neg <- s1s2_neg[!(duplicated(s1s2_neg[, "locus_tag_seed1"]) & duplicated(s1s2_neg[, "locus_tag_seed2"])), ]
x <- s1s2_neg[, "bestSNP_lod_seed1"]
y <- s1s2_neg[, "bestSNP_lod_seed2"]
ks.test(x, rnorm(length(x), mean(x), sd(x)))
ks.test(y, rnorm(length(y), mean(y), sd(y)))
shapiro.test(x)
shapiro.test(y)

cor_test <- cor.test(x, y, method = "spearman")
p_s1s2_neg <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05)
    scale_x_continuous(limits=c(5.0, 40)) + 
    scale_y_continuous(limits=c(5.0, 40)) + 
    geom_rug(col = rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("highest LOD, seed (repl. 1)") + ylab("highest LOD, seed (repl. 2)") +
    ggtitle(
        paste(expression(~ rho), "=", round(cor_test$estimate[["rho"]], 3))) +
    theme_bw()
ggsave(p, file = "plot_neg_seed1_seed2_scatter.pdf", device = "pdf", width = 16, height = 16.8, units = "cm", useDingbats=FALSE)

## seed1/leaf2
s1l2_neg <- map_neg[map_neg[, "locus_tag_seed1"] != "" & map_neg[, "locus_tag_leaf2"] != "", ]
##s1l2_neg <- s1l2_neg[!(duplicated(s1l2_neg[, "locus_tag_seed1"]) & duplicated(s1l2_neg[, "locus_tag_leaf2"])), ]
x <- s1l2_neg[, "bestSNP_lod_seed1"]
y <- s1l2_neg[, "bestSNP_lod_leaf2"]
ks.test(x, rnorm(length(x), mean(x), sd(x)))
ks.test(y, rnorm(length(y), mean(y), sd(y)))
shapiro.test(x)
shapiro.test(y)

cor_test <- cor.test(x, y, method = "spearman")
p <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits=c(5.0, 40)) + 
    scale_y_continuous(limits=c(5.0, 40)) + 
    geom_rug(col=rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("seed (replicate 1)") + ylab("leaf") +
    ggtitle(paste("LOD uniq loci neg", cor_test$estimate, cor_test$p.value)) + theme_bw()
ggsave(p, file = "plot_neg_seed1_leaf2_scatter.pdf", device = "pdf", width = 16, height = 16.8, units = "cm", useDingbats=FALSE)

## seed2/leaf2
s2l2_neg <- map_neg[map_neg[, "locus_tag_seed2"] != "" & map_neg[, "locus_tag_leaf2"] != "", ]
##s2l2_neg <- s2l2_neg[!(duplicated(s2l2_neg[, "locus_tag_seed1"]) & duplicated(s2l2_neg[, "locus_tag_seed2"])), ]
x <- s2l2_neg[, "bestSNP_lod_seed2"]
y <- s2l2_neg[, "bestSNP_lod_leaf2"]
ks.test(x, rnorm(length(x), mean(x), sd(x)))
ks.test(y, rnorm(length(y), mean(y), sd(y)))
shapiro.test(x)
shapiro.test(y)

cor_test <- cor.test(x, y, method = "spearman")
p <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits=c(5.0, 40)) + 
    scale_y_continuous(limits=c(5.0, 40)) + 
    geom_rug(col=rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("seed (replicate 2)") + ylab("leaf") +
    ggtitle(paste("LOD uniq loci neg", cor_test$estimate, cor_test$p.value)) + theme_bw()
ggsave(p, file = "plot_neg_seed2_leaf2_scatter.pdf", device = "pdf", width = 16, height = 16.8, units = "cm", useDingbats=FALSE)

## positive
## seed1/seed2
s1s2_pos <- map_pos[map_pos[, "locus_tag_seed1"] != "" & map_pos[, "locus_tag_seed2"] != "", ]
##s1s2_pos <- s1s2_pos[!(duplicated(s1s2_pos[, "locus_tag_seed1"]) & duplicated(s1s2_pos[, "locus_tag_seed2"])), ]
x <- s1s2_pos[, "bestSNP_lod_seed1"]
y <- s1s2_pos[, "bestSNP_lod_seed2"]
ks.test(x, rnorm(length(x), mean(x), sd(x)))
ks.test(y, rnorm(length(y), mean(y), sd(y)))
shapiro.test(x)
shapiro.test(y)

cor_test <- cor.test(x, y, method = "spearman")
p <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
  scale_x_continuous(limits=c(5.0, 40)) + 
  scale_y_continuous(limits=c(5.0, 40)) + 
  geom_rug(col=rgb(.5, 0, 0), alpha = 0.03) + 
  xlab("seed (replicate 1)") + ylab("seed (replicate 2)") +
  ggtitle(paste("LOD uniq loci pos", cor_test$estimate, cor_test$p.value)) + theme_bw()
ggsave(p, file = "plot_pos_seed1_seed2_scatter.pdf", device = "pdf", width = 16, height = 16.8, units = "cm", useDingbats=FALSE)

## seed1/leaf2
s1l2_pos <- map_pos[map_pos[, "locus_tag_seed1"] != "" & map_pos[, "locus_tag_leaf2"] != "", ]
##s1l2_pos <- s1l2_pos[!(duplicated(s1l2_pos[, "locus_tag_seed1"]) & duplicated(s1l2_pos[, "locus_tag_leaf2"])), ]
x <- s1l2_pos[, "bestSNP_lod_seed1"]
y <- s1l2_pos[, "bestSNP_lod_leaf2"]
ks.test(x, rnorm(length(x), mean(x), sd(x)))
ks.test(y, rnorm(length(y), mean(y), sd(y)))
shapiro.test(x)
shapiro.test(y)

cor_test <- cor.test(x, y, method = "spearman")
p <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
    scale_x_continuous(limits=c(5.0, 40)) + 
    scale_y_continuous(limits=c(5.0, 40)) + 
    geom_rug(col=rgb(.5, 0, 0), alpha = 0.03) + 
    xlab("seed (replicate 1)") + ylab("leaf") +
    ggtitle(paste("LOD uniq loci pos", cor_test$estimate, cor_test$p.value)) + theme_bw()
ggsave(p, file = "plot_pos_seed1_leaf2_scatter.pdf", device = "pdf", width = 16, height = 16.8, units = "cm", useDingbats=FALSE)

## seed2/leaf2
s2l2_pos <- map_pos[map_pos[, "locus_tag_seed2"] != "" & map_pos[, "locus_tag_leaf2"] != "", ]
##s2l2_pos <- s2l2_pos[!(duplicated(s2l2_pos[, "locus_tag_seed2"]) & duplicated(s2l2_pos[, "locus_tag_leaf2"])), ]
x <- s2l2_pos[, "bestSNP_lod_seed2"]
y <- s2l2_pos[, "bestSNP_lod_leaf2"]
ks.test(x, rnorm(length(x), mean(x), sd(x)))
ks.test(y, rnorm(length(y), mean(y), sd(y)))
shapiro.test(x)
shapiro.test(y)

cor_test <- cor.test(x, y, method = "spearman")
p <- qplot(x, y, data = data.frame(x = x, y = y), alpha = .05) +
  scale_x_continuous(limits=c(5.0, 40)) + 
  scale_y_continuous(limits=c(5.0, 40)) + 
  geom_rug(col=rgb(.5, 0, 0), alpha = 0.03) + 
  xlab("seed (replicate 2)") + ylab("leaf") +
  ggtitle(paste("LOD uniq loci pos", cor_test$estimate, cor_test$p.value)) + theme_bw()
ggsave(p, file = "plot_pos_seed2_leaf2_scatter.pdf", device = "pdf", width = 16, height = 16.8, units = "cm", useDingbats=FALSE)

###############################################################################
## UpSetR

map_neg <- read.table(
    "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_neg_cut.txt", 
    header = TRUE, sep = "\t")
map_pos <- read.table(
    "~/GitHub/GWAS_arabidopsis_seed/gwas_complete_met_all_trueLociLOD_pos_cut.txt", 
    header = TRUE, sep = "\t")



library(UpSetR)
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), header =T, sep=";")
mutations <- read.csv(system.file("extdata", "mutations.csv", package = "UpSetR"), header = T, sep = ",")
upset(movies)
upset(mutations)

## negative
up <- map_neg
up <- up[, c("locus_tag_seed1", "locus_tag_seed2", "locus_tag_leaf2", "met_rep1", "met_rep2")]
up[, 1] <- ifelse(up[, 1] == "", 0, 1)
up[, 2] <- ifelse(up[, 2] == "", 0, 1)
up[, 3] <- ifelse(up[, 3] == "", 0, 1)
pdf("plot_upset_neg.pdf")
upset(up[, 1:3])
dev.off()




## the trueLociLOD will now be truncated such that they only contain the mass
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
cols <- c("locus_tag_seed1", "locus_tag_seed2", "locus_tag_leaf2")
##cols_lod <- c("bestSNP_lod_seed1", "bestSNP_lod_seed2", "bestSNP_lod_leaf2")
binary_mat <- trueLociLOD |>
    dplyr::filter(met_rep1 %in% mapping_df$rep1) |>
    dplyr::filter(met_rep2 %in% mapping_df$rep2) |>
    dplyr::select(all_of(cols))
binary_mat <- ifelse(is.na(binary_mat) | binary_mat == "", 0, 1) |>
    as.data.frame()
l <- lapply(binary_mat, function(col_i) {
    .names <- ifelse(col_i == 1, seq_len(nrow(binary_mat)), 0)
    .names[.names != 0]
})
l <- UpSetR::fromList(l)
## for normalized
cols <- c("locus_tag_seed1", "locus_tag_seed2", "locus_tag_leaf2", 
          "locus_tag_leaf_feng")
binary_mat_normalized <- trueLociLOD_normalized |>
    dplyr::filter(
        met_rep1 %in% mapping_df_all$rep1 | is.na(met_rep1)) |>
    dplyr::filter(
        met_rep2 %in% mapping_df_all$rep2 | is.na(met_rep2)) |>
    dplyr::filter(
        met_leaf_feng %in% mapping_df_all$rep_feng | is.na(met_leaf_feng)) |>
    dplyr::select(all_of(cols))
binary_mat_normalized <- ifelse(
    is.na(binary_mat_normalized) | binary_mat_normalized == "", 0, 1) |>
    as.data.frame()
l_normalized <- lapply(binary_mat_normalized, function(col_i) {
    .names <- ifelse(col_i == 1, seq_len(nrow(binary_mat_normalized)), 0)
    .names[.names != 0]
})
l_normalized <- UpSetR::fromList(l_normalized)
## for batch
cols <- c("locus_tag_seed1", "locus_tag_seed2", "locus_tag_leaf2", 
          "locus_tag_leaf_feng")
binary_mat_batch <- trueLociLOD_batch |>
    dplyr::filter(
        met_rep1 %in% mapping_df_all$rep1 | is.na(met_rep1)) |>
    dplyr::filter(
        met_rep2 %in% mapping_df_all$rep2 | is.na(met_rep2)) |>
    dplyr::filter(
        met_leaf_feng %in% mapping_df_all$rep_feng | is.na(met_leaf_feng)) |>
    dplyr::select(all_of(cols))
binary_mat_batch <- ifelse(
    is.na(binary_mat_batch) | binary_mat_batch == "", 0, 1) |>
    as.data.frame()
l_batch <- lapply(binary_mat_batch, function(col_i) {
    .names <- ifelse(col_i == 1, seq_len(nrow(binary_mat_batch)), 0)
    .names[.names != 0]
})
l_batch <- UpSetR::fromList(l_batch)

## do the actual plotting
UpSetR::upset(l, order.by = "freq")
UpSetR::upset(l_normalized, order.by = "freq")
UpSetR::upset(l_batch, order.by = "freq")



















##################################################################################

## only leaf2
df <- data.frame(values = table(up[up[, 1] == 0 & up[, 2] == 0 & up[, 3] == 1, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_neg_number_mass_features_l2.pdf", device = "pdf", useDingbats=FALSE)

## only seed1
df <- data.frame(values = table(up[up[, 1] == 1 & up[, 2] == 0 & up[, 3] == 0, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_neg_number_mass_features_s1.pdf", device = "pdf", useDingbats=FALSE)

## only seed2
df <- data.frame(values = table(up[up[, 1] == 0 & up[, 2] == 1 & up[, 3] == 0, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_neg_number_mass_features_s2.pdf", device = "pdf", useDingbats=FALSE)

## seed1&seed2
df <- data.frame(values = table(up[up[, 1] == 1 & up[, 2] == 1 & up[, 3] == 0, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_neg_number_mass_features_s1s2.pdf", device = "pdf", useDingbats=FALSE)

## seed1&leaf2
df <- data.frame(values = table(up[up[, 1] == 1 & up[, 2] == 0 & up[, 3] == 1, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_neg_number_mass_features_s1l2.pdf", device = "pdf", useDingbats=FALSE)

## seed2&leaf2
df <- data.frame(values = table(up[up[, 1] == 0 & up[, 2] == 1 & up[, 3] == 1, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_neg_number_mass_features_s2l2.pdf", device = "pdf", useDingbats=FALSE)

## seed1&seed2&leaf2
df <- data.frame(values = table(up[up[, 1] == 1 & up[, 2] == 1 & up[, 3] == 1, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_neg_number_mass_features_s1s2l2.pdf", device = "pdf", useDingbats=FALSE)

## positive
up <- map_pos
up <- up[, c("locus_tag_seed1", "locus_tag_seed2", "locus_tag_leaf2", "met_rep1", "met_rep2")]
up[, 1] <- ifelse(up[, 1] == "", 0, 1)
up[, 2] <- ifelse(up[, 2] == "", 0, 1)
up[, 3] <- ifelse(up[, 3] == "", 0, 1)
pdf("plot_upset_pos.pdf")
upset(up[, 1:3])
dev.off()


## only leaf2
df <- data.frame(values = table(up[up[, 1] == 0 & up[, 2] == 0 & up[, 3] == 1, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_pos_number_mass_features_l2.pdf", device = "pdf", useDingbats=FALSE)

## only seed1
df <- data.frame(values = table(up[up[, 1] == 1 & up[, 2] == 0 & up[, 3] == 0, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_pos_number_mass_features_s1.pdf", device = "pdf", useDingbats=FALSE)

## only seed2
df <- data.frame(values = table(up[up[, 1] == 0 & up[, 2] == 1 & up[, 3] == 0, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_pos_number_mass_features_s2.pdf", device = "pdf", useDingbats=FALSE)

## seed1&seed2
df <- data.frame(values = table(up[up[, 1] == 1 & up[, 2] == 1 & up[, 3] == 0, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_pos_number_mass_features_s1s2.pdf", device = "pdf", useDingbats=FALSE)

## seed1&leaf2
df <- data.frame(values = table(up[up[, 1] == 1 & up[, 2] == 0 & up[, 3] == 1, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_pos_number_mass_features_s1l2.pdf", device = "pdf", useDingbats=FALSE)

## seed2&leaf2
df <- data.frame(values = table(up[up[, 1] == 0 & up[, 2] == 1 & up[, 3] == 1, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_pos_number_mass_features_s2l2.pdf", device = "pdf", useDingbats=FALSE)

## seed1&seed2&leaf2
df <- data.frame(values = table(up[up[, 1] == 1 & up[, 2] == 1 & up[, 3] == 1, "met_rep1"]))
p <- ggplot(df, aes(x = values.Freq)) + geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(c(0,1)) +xlim(c(0,75)) + theme_bw()
ggsave(p, file = "plot_distribution_pos_number_mass_features_s1s2l2.pdf", device = "pdf", useDingbats=FALSE)

