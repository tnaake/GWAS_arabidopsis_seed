
## add columns for the file objects that do not have all columns
add_column <- function(file_n, cols = cols) {
    
    ## get cols that are not present in file_n
    .cols <- cols[!cols %in% colnames(file_n)]
    
    ## create data.frame and fill with "." for columns that are not present
    ## in file_n
    .df <- matrix(".", nrow = nrow(file_n), ncol = length(.cols)) |>
        as.data.frame()
    names(.df) <- .cols
    
    ## add the data.frame to .df
    file_n <- cbind(file_n, .df)  
    
    ## sort the file_n object according to cols and return
    file_n[, cols]
}

combine_files <- function(path, file = "Thomasgwas_gene_info_final.txt") {
    f <- list.files(path = path, pattern = "gene_info.txt", recursive = FALSE,
        full.names = TRUE)
    f <- lapply(f, function(f_i) read.delim(f_i, sep = "\t", dec = ".", 
        quote = ""))
    
    ## get unique colnames
    cols <- lapply(f, function(f_i) colnames(f_i)) |>
        unlist() |>
        unique()
    
    ## add columns for the file objects that do not have all columns
    f_added <- lapply(f, function(f_i) add_column(file_n = f_i, cols = cols))
    
    ## combine the file objects
    file_final <- do.call("rbind", f_added)
    
    write.table(file_final, file = file, 
        quote = FALSE, sep = "\t", dec = ".", row.names = FALSE)
    
}


## old files
setwd("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/files_Feng/old")
path <- getwd()
file <- "Thomasgwas_gene_info_final.txt"
combine_files(path = path, file = file)

## normalized
setwd("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/files_Feng/normalized")
path <- getwd()
file <- "Thomasgwas_gene_info_final_normalized.txt"
combine_files(path = path, file = file)

## normalized_batch_corrected
setwd("~/GitHub/GWAS_arabidopsis_seed/GWAS_leaf_rep2/files_Feng/normalized_batch_corrected/")
path <- getwd()
file <- "Thomasgwas_gene_info_final_normalized_batch_corrected.txt"
combine_files(path = path, file = file)
