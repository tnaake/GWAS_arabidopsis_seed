## replicate 1, negative
## set working directory

classes_samples <- c(rep("C24", 4), rep("Col_0", 10), rep("SALK_008908_3", 20),
    rep("SALK_011180_5", 5), rep("SALK_020876_4", 5), rep("SALK_021216_4", 5),
    rep("SALK_024438_5", 5), rep("SALK_027837_5", 5), rep("SALK_037430_4", 5),
    rep("SALK_049338_4", 5), rep("SALK_072964_5", 5), rep("SALK_081021_6", 20),
    rep("SALK_201809C_2", 5), rep("SALK_203337C_2", 5), 
    rep("SALK_203919C_2", 20), rep("SALK_204674C_3", 5)
)


setwd("~/GitHub/GWAS_arabidopsis_seed/TDNA_rep1/neg/")

## load xcms
library(xcms)
xset_neg_rep1 <- xcmsSet(file = "./", method = "centWave", ppm = 5,
    snthresh = 10, peakwidth = c(5, 20), prefilter = c(1, 5000))
sampclass(xset_neg_rep1) <- classes_samples
xset2_neg_rep1 <- group(xset_neg_rep1, method = "density", minfrac = 0.5, 
    minsamp = 1, bw = 2, mzwid = 0.015)
xset3_neg_rep1 <- retcor(xset2_neg_rep1, family = "s", plottype = "m", 
    missing = 1, extra = 1, span = 1)
xset4_neg_rep1 <- group(xset3_neg_rep1, method = "density", minfrac = 0.5, 
    minsamp = 1, bw = 2, mzwid = 0.015)
xset5_neg_rep1 <- fillPeaks(xset4_neg_rep1, method = "chrom")
save("xset_neg_rep1", "xset2_neg_rep1", "xset3_neg_rep1", "xset4_neg_rep1", 
    "xset5_neg_rep1", 
    file = "./TDNA_replicate1_neg_xcms.RData")

## CAMERA
library(CAMERA)
an_neg_rep1 <- xsAnnotate(xset5_neg_rep1)
anF_neg_rep1 <- groupFWHM(an_neg_rep1, perfwhm = 0.6)
anI_neg_rep1 <- findIsotopes(anF_neg_rep1, mzabs = 0.01)
anIC_neg_rep1 <- groupCorr(anI_neg_rep1, cor_eic_th = 0.75, graphMethod = "lpc")
anFA_neg_rep1 <- findAdducts(anIC_neg_rep1, polarity = "negative")
pl_neg_rep1 <- getPeaklist(anFA_neg_rep1)
save("an_neg_rep1", "anF_neg_rep1", "anI_neg_rep1", "anIC_neg_rep1", 
    "anFA_neg_rep1", "pl_neg_rep1", 
    file = "./TDNA_replicate1_neg_CAMERA.RData")

write.table(pl_neg_rep1, file = "TDNA_replicate1_neg_peaklist.txt", sep = "\t",
    dec = ".", row.names = FALSE, quote = FALSE)

## replicate 1, positive mode
## set working directory
setwd("~/GitHub/GWAS_arabidopsis_seed/TDNA_rep1/pos/")

## load xcms
library(xcms)
xset_pos_rep1 <- xcmsSet(file = "./", method = "centWave", ppm = 5,
    snthresh = 10, peakwidth = c(5, 20), prefilter = c(1, 5000))
sampclass(xset_pos_rep1) <- classes_samples
xset2_pos_rep1 <- group(xset_pos_rep1, method = "density", minfrac = 0.5, 
    minsamp = 1, bw = 2, mzwid = 0.015)
xset3_pos_rep1 <- retcor(xset2_pos_rep1, family = "s", plottype = "m", 
    missing = 1, extra = 1, span = 1)
xset4_pos_rep1 <- group(xset3_pos_rep1, method = "density", minfrac = 0.5, 
    minsamp = 1, bw = 2, mzwid = 0.015)
xset5_pos_rep1 <- fillPeaks(xset4_pos_rep1, method = "chrom")
save("xset_pos_rep1", "xset2_pos_rep1", "xset3_pos_rep1", "xset4_pos_rep1", 
    "xset5_pos_rep1", 
    file = "./TDNA_replicate1_pos_xcms.RData")

## CAMERA
library(CAMERA)
an_pos_rep1 <- xsAnnotate(xset5_pos_rep1)
anF_pos_rep1 <- groupFWHM(an_pos_rep1, perfwhm = 0.6)
anI_pos_rep1 <- findIsotopes(anF_pos_rep1, mzabs = 0.01)
anIC_pos_rep1 <- groupCorr(anI_pos_rep1, cor_eic_th = 0.75, graphMethod = "lpc")
anFA_pos_rep1 <- findAdducts(anIC_pos_rep1, polarity = "positive")
pl_pos_rep1 <- getPeaklist(anFA_pos_rep1)
save("an_pos_rep1", "anF_pos_rep1", "anI_pos_rep1", "anIC_pos_rep1", 
     "anFA_pos_rep1", "pl_pos_rep1", 
     file = "./TDNA_replicate1_pos_CAMERA.RData")

write.table(pl_pos_rep1, file = "TDNA_replicate1_pos_peaklist.txt", sep = "\t",
            dec = ".", row.names = FALSE, quote = FALSE)


#####################
## replicate 2, negative mode
## set working directory
setwd("~/GitHub/GWAS_arabidopsis_seed/TDNA_rep2/neg/")

## load xcms
library(xcms)
xset_neg_rep2 <- xcmsSet(file = "./", method = "centWave", ppm = 5,
    snthresh = 10, peakwidth = c(5, 20), prefilter = c(1, 5000))
sampclass(xset_neg_rep2) <- classes_samples
xset2_neg_rep2 <- group(xset_neg_rep2, method = "density", minfrac = 0.5, 
    minsamp = 1, bw = 2, mzwid = 0.015)
xset3_neg_rep2 <- retcor(xset2_neg_rep2, family = "s", plottype = "m", 
    missing = 1, extra = 1, span = 1)
xset4_neg_rep2 <- group(xset3_neg_rep2, method = "density", minfrac = 0.5, 
    minsamp = 1, bw = 2, mzwid = 0.015)
xset5_neg_rep2 <- fillPeaks(xset4_neg_rep2, method = "chrom")
save("xset_neg_rep2", "xset2_neg_rep2", "xset3_neg_rep2", "xset4_neg_rep2", 
     "xset5_neg_rep2", 
     file = "./TDNA_replicate2_neg_xcms.RData")

## CAMERA
library(CAMERA)
an_neg_rep2 <- xsAnnotate(xset5_neg_rep2)
anF_neg_rep2 <- groupFWHM(an_neg_rep2, perfwhm = 0.6)
anI_neg_rep2 <- findIsotopes(anF_neg_rep2, mzabs = 0.01)
anIC_neg_rep2 <- groupCorr(anI_neg_rep2, cor_eic_th = 0.75, graphMethod = "lpc")
anFA_neg_rep2 <- findAdducts(anIC_neg_rep2, polarity = "negative")
pl_neg_rep2 <- getPeaklist(anFA_neg_rep2)
save("an_neg_rep2", "anF_neg_rep2", "anI_neg_rep2", "anIC_neg_rep2", 
    "anFA_neg_rep2", "pl_neg_rep2", 
    file = "./TDNA_replicate2_neg_CAMERA.RData")

write.table(pl_neg_rep2, file = "TDNA_replicate2_neg_peaklist.txt", sep = "\t",
            dec = ".", row.names = FALSE, quote = FALSE)
