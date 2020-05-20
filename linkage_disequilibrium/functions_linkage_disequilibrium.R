# arrange x position for gene labels in a manhattan plot
arrange <- function(x, diff = 1.5e7, rg){
	  
    n <- length(x)
	  
	  if(n == 1)
		    return(x)
	  
	  if(n %% 2 == 0) {
	    
		    if(x[n / 2 + 1] - x[n / 2] < diff) {
			      d <- (diff - (x[n / 2 + 1] - x[n / 2])) / 2
			      x[n / 2 + 1] <- x[n / 2 + 1] + d
			      x[n / 2] <- x[n / 2] - d
		    }
	  }
	  k <- if(n %% 2 == 0) n / 2 + 2 else (n + 1) / 2 + 1
	  while(k <= n) {
		    if(x[k] - x[k - 1] < diff)
			      x[k] <- x[k - 1] + diff
		    k <- k + 1
	  }
	  k <- if(n %% 2 == 0) n / 2 - 1 else (n + 1) / 2 - 1
	  while (k > 0) {
		    if(x[k + 1] - x[k] < diff)
			      x[k] <- x[k + 1] - diff
		    k <- k - 1
	  }
	  if(x[1] - rg[1] < diff / 2) {
		    x[1] <- diff / 2
		    for(k in 2:n) {
			      if(x[k] - x[k - 1] < diff)
				        x[k] <- x[k - 1] + diff
		    }
	  }
	  if(x[n] > rg[2] - diff / 2) {
		    x[n] <- rg[2] - diff / 2
		    for(k in (n - 1):1) {
			      if(x[k + 1] - x[k] < diff)
				        x[k] <- x[k + 1] - diff
		    }
	  }
	  return(x)
}

# make plot
makeLDplot <- function(gwas, tair, snp, genes, outputFile, plotTitle, lodFilter, lodMax, winUp, winDown) {    

    gwas$lod <- -log10(gwas$P.value)
    gwas <- gwas[gwas$lod > lodFilter, ]
    gwas <- gwas[order(gwas$lod), ] ## put less significant SNP first

    lod <- gwas$lod
    nfo <- tair[match(names(genes), tair$locus_tag), c("chrom", "pos_i", "pos_f")]

    tmp <-  which(tair$pos_i >= min(nfo$pos_i) & tair$pos_f <= max(nfo$pos_f)
        & tair$chrom == nfo[1, 1])

    nfo2 <- tair[tmp, c("chrom", "pos_i", "pos_f", "strand", "locus_tag")]
    nfo2 <- nfo2[!nfo2$locus_tag %in% names(genes), ]

    j <- which(gwas$Chromosome == nfo[1, 1] & gwas$Pos >= min(nfo$pos_i) - winDown &
            gwas$Pos <= max(nfo$pos_f) + winUp)

    z <- cbind(gwas$Pos[j] / 1e6, lod[j])

    jj <- which.max(lod[j])
    stopifnot(jj == nrow(z))
    ld.mat <- r2fast(snp, gwas$SNP[j])
    ld.mat <- ld.mat[, jj]
    ld.mat[is.na(ld.mat)] <- 1

    ## define colors
    col <- brewer.pal(10, "Spectral")
    col[1] <- "white"
    ld.col <- col[cut(ld.mat, c(-1, 1:9 / 10, 1.1), labels = FALSE)]

    ## axis limits
    ylab <- expression(- "log"[10] (P))
    ymax <- lodMax
    xlim <- c( min(nfo$pos_i) - winDown ,  max(nfo$pos_f) + winUp) / 1e6

    ## actual plot
    pdf(file = outputFile)
        
    m <- matrix(c(1, 1, 2, 0), 2, 2)
    layout(m, width=c(1, lcm(2)), height=c(lcm(8), 1))
    par(mar=c(5.1, 5, 5, 1)) 

    panel.first <- function() {
        ## plot horizontal line (currently commented)
        ## abline(v=unlist(nfo[,2:3])/1e6, col='gray')
        rect(nfo[, 2] / 1e6, -5, nfo[, 3] / 1e6, lodMax + 2.5, col = 'gray90')
        abline(h = 5.3, col = 'gray', lty = 2)
    }

    plot(z[-jj,], pch = 21, las = 1, ylab = ylab, xlim = xlim, ylim = c(0, ymax), 
        bg = ld.col[-jj], cex = 1.5,  main = plotTitle,
        xlab = sprintf("Chr. %d (MB)", nfo[1, 1]), panel.first = panel.first())
    ## plot SNP with highest LOD in red and as a diamond
    points(z[jj, 1], z[jj, 2], pch = 18, col = 'red', cex = 2) 

    # plot the R^2 color key
    par(mar=c(0, 0, 5, 3)) 
    image(1, 1:10, matrix(1:10, nrow = 1), xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "", col = col)
    axis(4, at = seq(0.5, 10.5, by = 2), 
        labels = format(seq(0, 1, length = 6)), las = 1)
    mtext(expression(r^2), 3, line = 1) 
    box()
    dev.off()
    invisible()
}
