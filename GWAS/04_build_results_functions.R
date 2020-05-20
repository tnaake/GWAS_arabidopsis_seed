## helper functions for 04_build_results.R 
overlap <- function(l1, l2, l3, l1_info = seed1, l2_info = seed2, l3_info = leaf2) {
    
    res <- matrix(0, nrow = 0, ncol = 9)
    
    for (k in seq_along(l1)) {
        l1_k <- l1[[k]]
        l1_k_b <- l1_k[1]
        l1_k_e <- l1_k[length(l1_k)]
        
        locusID_l1 <- l1_info[l1_info[, "locus_tag"] == l1_k_b, "locusID"]
        bestSNP_lod_l1 <- l1_info[l1_info[, "locus_tag"] == l1_k_b, "best_SNP_lod"]
        
        ## check all combinations of l2 --> l3
        for (m in seq_along(l2)) {
            l2_m <- l2[[m]]
            l2_m_b <- l2_m[1]
            l2_m_e <- l2_m[length(l2_m)]
            
            locusID_l2 <- l2_info[l2_info[, "locus_tag"] == l2_m_b, "locusID"]
            bestSNP_lod_l2 <- l2_info[l2_info[, "locus_tag"] == l2_m_b, "best_SNP_lod"]
            
            if (any(l1_k %in% l2_m)) { 
                
                for (n in seq_along(l3)) {
                    l3_n <- l3[[n]]
                    l3_n_b <- l3_n[1]
                    l3_n_e <- l3_n[length(l3_n)]
                    
                    locusID_l3 <- l3_info[l3_info[, "locus_tag"] == l3_n_b, "locusID"]
                    bestSNP_lod_l3 <- l3_info[l3_info[, "locus_tag"] == l3_n_b, "best_SNP_lod"]
                    
                    if (any(l1_k %in% l3_n)) { ## overlap between l1, l2, l3
                        res <- rbind(res, c(locusID_l1, bestSNP_lod_l1, 
                            paste(l1_k_b, l1_k_e, sep = ".."), 
                            locusID_l2, bestSNP_lod_l2, 
                            paste(l2_m_b, l2_m_e, sep = ".."), 
                            locusID_l3, bestSNP_lod_l3, 
                            paste(l3_n_b, l3_n_e, sep = "..")))
                    } else { ## only overlap between l1 and l2 
                        res <- rbind(res, c(locusID_l1, bestSNP_lod_l1, 
                            paste(l1_k_b, l1_k_e, sep = ".."), 
                            locusID_l2, bestSNP_lod_l2, 
                            paste(l2_m_b, l2_m_e, sep = ".."), "", "", ""))
                    }
                }
                if (all(unlist(l3) == ""))
                    res <- rbind(res, c(locusID_l1, bestSNP_lod_l1, 
                        paste(l1_k_b, l1_k_e, sep = ".."), locusID_l2, 
                        bestSNP_lod_l2, paste(l2_m_b, l2_m_e, sep = ".."), "", "", ""))
            } else {
                res <- rbind(res, c(locusID_l1, bestSNP_lod_l1, 
                    paste(l1_k_b, l1_k_e, sep = ".."), "", "", "", "", "", ""))    
            }
        }
        if (all(unlist(l2) == ""))
            res <- rbind(res, c(locusID_l1, bestSNP_lod_l1, 
                paste(l1_k_b, l1_k_e, sep = ".."), "", "", "", "", "", ""))
        
        ## check all combinations of l3 --> l2 
        for (m in seq_along(l3)) {
            l3_m <- l3[[m]]
            l3_m_b <- l3_m[1]
            l3_m_e <- l3_m[length(l3_m)]
            
            locusID_l3 <- l3_info[l3_info[, "locus_tag"] == l3_m_b, "locusID"]
            bestSNP_lod_l3 <- l3_info[l3_info[, "locus_tag"] == l3_m_b, "best_SNP_lod"]
            
            if (any(l1_k %in% l3_m)) { 
                
                for (n in seq_along(l2)) {
                    l2_n <- l2[[n]]
                    l2_n_b <- l2_n[1]
                    l2_n_e <- l2_n[length(l2_n)]
                    
                    locusID_l2 <- l2_info[l2_info[, "locus_tag"] == l2_n_b, "locusID"]
                    bestSNP_lod_l2 <- l2_info[l2_info[, "locus_tag"] == l2_n_b, "best_SNP_lod"]
                    
                    if (any(l1_k %in% l2_n)) { ## overlap between l1, l2, l3
                        res <- rbind(res, c(locusID_l1, bestSNP_lod_l1, 
                            paste(l1_k_b, l1_k_e, sep = ".."), locusID_l2, 
                            bestSNP_lod_l2, paste(l2_n_b, l2_n_e, sep = ".."), 
                            locusID_l3, bestSNP_lod_l3, paste(l3_m_b, l3_m_e, sep = "..")))
                    } else { ## only overlap between l1 and l3
                        res <- rbind(res, c(locusID_l1, bestSNP_lod_l1, 
                            paste(l1_k_b, l1_k_e, sep = ".."), "", "", "", 
                            locusID_l3, bestSNP_lod_l3, paste(l3_m_b, l3_m_e, sep = "..")))
                    } 
                }
                if (all(unlist(l2) == ""))
                    res <- rbind(res, c(locusID_l1, bestSNP_lod_l1, 
                        paste(l1_k_b, l1_k_e, sep = ".."), "", "", "", 
                        locusID_l3, bestSNP_lod_l3, paste(l3_m_b, l3_m_e, sep = "..")))
            } else {
                res <- rbind(res, c(locusID_l1, bestSNP_lod_l1, 
                    paste(l1_k_b, l1_k_e, sep = ".."), "", "", "", "", "", ""))
            }
        }
        if (all(unlist(l3) == ""))
            res <- rbind(res, c(locusID_l1, bestSNP_lod_l1, 
                paste(l1_k_b, l1_k_e, sep = ".."), "", "", "", "", "", ""))
    }
    
    ## remove duplicates
    dup_l1 <- duplicated(res[, 3])
    dup_l2 <- duplicated(res[, 6])
    dup_l3 <- duplicated(res[, 9])
    
    ## remove duplicates (all identical)
    res <- res[!(dup_l1 & dup_l2 & dup_l3), ]
    res <- matrix(res, ncol = 9)
    
    ## remove those that have one empty entry but are other-wise identical
    inds <- which(res == "", arr.ind = TRUE)
    rem <- logical(nrow(res))
    for (j in seq_along(inds[, "row"])) {
        res_log <- res[, -inds[j, "col"]] == res[inds[j, "row"], -inds[j, "col"]]
        if (is.matrix(res_log)) {
            res_log_all <- apply(res_log, 1, function(x) all(x))
            ## check now TRUE rows --> remove the row that has one "" at the 
            ## position inds[j, "col"]
            if (sum(res_log_all) > 1) 
                rem[res[, inds[j, "col"]] == "" & res_log_all] <- TRUE
        }
    }
    
    res <- res[!rem, ]
    
    res <- matrix(res, ncol = 9)
    return(res)   
}

overlap_l <- function(l) {
    l1 <- l$l1
    l2 <- l$l2
    l3 <- l$l3
    
    if (all(unlist(l1) != "")) {
        ol_1 <- overlap(l1, l2, l3, l1_info = seed1, l2_info = seed2, l3_info = leaf2)    
    } else {
        ol_1 <- matrix("", ncol = 9, nrow = 0)
    }
    colnames(ol_1) <- c("locusID_seed1", "bestSNP_lod_seed1", "locus_tag_seed1", 
        "locusID_seed2", "bestSNP_lod_seed2", "locus_tag_seed2", 
        "locusID_leaf2", "bestSNP_lod_leaf2", "locus_tag_leaf2")
    
    if (all(unlist(l2) != "")) {
        ol_2 <- overlap(l2, l1, l3, l1_info = seed2, l2_info = seed1, l3_info = leaf2)    
    } else {
        ol_2 <- matrix("", ncol = 9, nrow = 0)
    }
    colnames(ol_2) <- c("locusID_seed2", "bestSNP_lod_seed2", "locus_tag_seed2", 
        "locusID_seed1", "bestSNP_lod_seed1", "locus_tag_seed1", 
        "locusID_leaf2", "bestSNP_lod_leaf2", "locus_tag_leaf2")
    
    if (all(unlist(l3) != "")) {
        ol_3 <- overlap(l3, l1, l2, l1_info = leaf2, l2_info = seed1, l3_info = seed2)    
    } else {
        ol_3 <- matrix("", ncol = 9, nrow = 0)
    }
    colnames(ol_3) <- c("locusID_leaf2", "bestSNP_lod_leaf2", "locus_tag_leaf2", 
        "locusID_seed1", "bestSNP_lod_seed1", "locus_tag_seed1", 
        "locusID_seed2", "bestSNP_lod_seed2", "locus_tag_seed2")
    
    ## rbind the results
    c_ord <- c("locusID_seed1", "bestSNP_lod_seed1", "locus_tag_seed1", 
        "locusID_seed2", "bestSNP_lod_seed2", "locus_tag_seed2", 
        "locusID_leaf2", "bestSNP_lod_leaf2", "locus_tag_leaf2")
    res <- rbind(ol_1[, c_ord], ol_2[, c_ord], ol_3[, c_ord])
    
    ## remove duplicates
    dup_l1 <- duplicated(res[, "locus_tag_seed1"])
    dup_l2 <- duplicated(res[, "locus_tag_seed2"])
    dup_l3 <- duplicated(res[, "locus_tag_leaf2"])
    
    ## remove duplicates (all identical)
    res <- res[!(dup_l1 & dup_l2 & dup_l3), ]
    res <- matrix(res, ncol = 9)
    
    ## remove those that have one empty entry but are other-wise identical
    inds <- which(res == "", arr.ind = TRUE)
    rem <- logical(nrow(res))
    for (j in seq_along(inds[, "row"])) {
        res_log <- res[, -inds[j, "col"]] == res[inds[j, "row"], -inds[j, "col"]]
        if (is.matrix(res_log)) {
            res_log_all <- apply(res_log, 1, function(x) all(x))
            
            ## check now TRUE rows --> remove the row that has one "" at the 
            ## position inds[j, "col"]
            if (sum(res_log_all) > 1) 
                rem[res[, inds[j, "col"]] == "" & res_log_all] <- TRUE    
        }
    }
    
    ## remove those that have two empty entries and combined are identical 
    ## to a two entried row
    row_table <- table(inds[, "row"])
    row_table <- as.numeric(names(which(row_table == 2)))
    
    if (length(row_table) > 0) {
        row_combn <- combn(row_table, 2)
        for (j in seq_len(ncol(row_combn))) {
            comb_row <- apply(res[row_combn[, j], ], 2, paste, collapse = "")
            res_comb <- t(apply(res, 1, function(x) x == comb_row))
            if (any(apply(res_comb, 1, function(x) all(x)))) {
                for (k in row_combn[, j])
                    res[k, ] <- comb_row
            }
        }
        
        dup_l1 <- duplicated(res[, "locus_tag_seed1"])
        dup_l2 <- duplicated(res[, "locus_tag_seed2"])
        dup_l3 <- duplicated(res[, "locus_tag_leaf2"])
        
        ## remove duplicates (all identical)
        rem[(dup_l1 & dup_l2 & dup_l3)] <- TRUE
    }
    ## finally remove and return
    res <- res[!rem, ]
    res <- matrix(res, ncol = 9)
    colnames(res) <- c_ord
    
    rem <- logical(nrow(res))
    
    ## remove those that have two empty entries when there is one row with two 
    ## entries
    dup_l1 <- duplicated(res[, "locus_tag_seed1"])
    dup_l2 <- duplicated(res[, "locus_tag_seed2"])
    dup_l3 <- duplicated(res[, "locus_tag_leaf2"])
    
    inds <- which(dup_l1 & res[, "locus_tag_seed1"] != "" | dup_l2 & res[, "locus_tag_seed2"] != "" | dup_l3 & res[, 3] != "")
    for (j in inds) {
        if (dup_l1[j] & res[j, "locus_tag_seed1"] != "") {
            tmp <- res[j, "locus_tag_seed1"]
            inds_j <- which(res[, "locus_tag_seed1"] == tmp)
            rem[inds_j[which.max(apply(res[inds_j, ] == "", 1, sum))]] <- TRUE
        }
        if (dup_l2[j] & res[j, "locus_tag_seed2"] != "") {
            tmp <- res[j, "locus_tag_seed2"]
            inds_j <- which(res[, "locus_tag_seed2"] == tmp)
            rem[inds_j[which.max(apply(res[inds_j, ] == "", 1, sum))]] <- TRUE
        }
        if (dup_l3[j] & res[j, "locus_tag_leaf2"] != "") {
            tmp <- res[j, "locus_tag_leaf2"]
            inds_j <- which(res[, "locus_tag_leaf2"] == tmp)
            rem[inds_j[which.max(apply(res[inds_j, ] == "", 1, sum))]] <- TRUE
        }
    }
    res <- res[!rem, ]

    return(res)
}

test_l1 <- list(l1 = list(c("a", "b", "c"), c("d", "e"), c("x", "y", "z")),
                l2 = list(c("b", "c"), c("d", "e"), c("f", "g")),
                l3 = list(c("d", "e"), c("f", "g")))

test_l2 <- list(l1 = list(c("a", "b", "c"), c("d", "e")),
                l2 = list(c("b", "c"), c("d", "e"), c("f", "g")),
                l3 = list(""))

test_l3 <- list(l1 = list(c("a", "b", "c"), c("d", "e")),
                l2 = list(""),
                l3 = list(c("x", "y", "z")))

test_l4 <- list(l1 = list(c("a", "b", "c"), c("d", "e")),
                l2 = list(),
                l3 = list())

test_l5 <- list(l1 = list(c("a", "b", "c"), c("d", "e")),
                l2 = list(c("f", "g"), c("h", "i")),
                l3 = list(c("j", "k"), c("l", "m")))

overlap_l(test_l1)
overlap_l(test_l2)
overlap_l(test_l3)
overlap_l(test_l4)
overlap_l(test_l5)
