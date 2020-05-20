setwd("~/winhome/Documents/01_GWAS/01_Data/GWAS_script/LD")

## script to make the LD figure
library(GenABEL)
library(RColorBrewer)
library(gplots)
load("snp_genabel.RData")
load("Data/TAIR9.RData")
source("functions_linkage_disequilibrium.R")
options(stringsAsFactors = FALSE)

## enter here the gene of interest
gene <- "AT5G07990" 

## enter here the metabolite of interest
cluster <- "Cluster_146144"

## define the genomic range (+- limit of genes around gene)
limit <- 50

## load GAPIT result file
gwas <- read.csv(paste("GAPIT.MLM.", cluster, ".GWAS.Results.csv", sep = ""))

## put here the genes you want to show (AGI code = gene name) / put the genes here that mapped to the feature
genes <- read.table("gene_list_tair9.txt", header = TRUE, sep = "\t", 
    quote = "\"", dec = ".", comment.char = "")

goi <- which(genes[, "locus_tag"] == gene) 

genes <- genes[(goi - limit):(goi + limit), 1]
names(genes) <- genes

## options for the plot (ranges)
lodFilter  <- 0
lodMax     <- 12.5 
name <- paste(cluster, gene, limit, sep = "_")
outputFile <- paste(name, "pdf", sep = ".")
plotTitle  <- name

## base pairs up and down stream of the genes of interest
winUp <- 5000
winDown <- 5000

## create the actual plot
makeLDplot(gwas, tair, snp, genes, outputFile, plotTitle, lodFilter, lodMax, winUp, winDown)
