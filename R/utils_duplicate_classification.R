

#' Classify gene pairs derived from segmental duplications
#'
#' @param anchor_pairs A 2-column data frame with anchor pairs in columns 1
#' and 2.
#' @param pairs A 2-column data frame with all duplicate pairs. This
#' is equivalent to the first 2 columns of the tabular output of BLAST-like
#' programs.
#'
#' @return A 3-column data frame with the variables:
#' \describe{
#'   \item{dup1}{Character, duplicated gene 1}
#'   \item{dup2}{Character, duplicated gene 2}
#'   \item{type}{Factor indicating duplication types, with levels
#'               "SD" (segmental duplication) or 
#'               "DD" (dispersed duplication).}
#' }
#' @rdname get_segmental
#' @export
#' @examples
#' data(diamond_intra)
#' data(yeast_annot)
#' data(yeast_seq)
#' blast_list <- diamond_intra
#' 
#' # Get processed annotation for S. cerevisiae
#' annotation <- syntenet::process_input(yeast_seq, yeast_annot)$annotation[1]
#' 
#' # Get list of intraspecies anchor pairs
#' anchor_pairs <- get_anchors_list(blast_list, annotation)
#' anchor_pairs <- anchor_pairs[[1]][, c(1, 2)]
#' 
#' # Get duplicate pairs from DIAMOND output
#' duplicates <- diamond_intra[[1]][, c(1, 2)]
#' dups <- get_segmental(anchor_pairs, duplicates)
get_segmental <- function(anchor_pairs = NULL, pairs = NULL) {
    
    names(pairs) <- c("dup1", "dup2")
    p <- pairs[pairs$dup1 != pairs$dup2, ]
    anchorp <- anchor_pairs
    
    duplicates <- p
    duplicates$type <- "DD"
    
    if(!is.null(anchorp)) {
        names(anchorp) <- c("anchor1", "anchor2")
        
        # Look for anchor pairs in duplicate pairs - vector-based approach
        p_vector <- paste0(p$dup1, p$dup2)
        anchor_vector <- c(
            paste0(anchorp$anchor1, anchorp$anchor2),
            paste0(anchorp$anchor2, anchorp$anchor1)
        )
        
        # Add column `type` with classification
        dup_mode <- ifelse(p_vector %in% anchor_vector, "SD", "DD")
        duplicates$type <- dup_mode
    }
    duplicates$type <- factor(duplicates$type, levels = c("SD", "DD"))
    
    return(duplicates)
}


#' Classify gene pairs derived from tandem and proximal duplications
#' 
#' @param pairs A 3-column data frame with columns \strong{dup1}, \strong{dup2},
#' and \strong{type} indicating duplicated gene 1, duplicated gene 2, and
#' the mode of duplication associated with the pair. This data frame
#' is returned by \code{get_segmental()}.
#' @param annotation_granges A processed GRanges object as in each element 
#' of the list returned by \code{syntenet::process_input()}.
#' @param proximal_max Numeric scalar with the maximum distance (in number
#' of genes) between two genes to consider them as proximal duplicates.
#' Default: 10.
#'
#' @return A 3-column data frame with the variables:
#' \describe{
#'   \item{dup1}{Character, duplicated gene 1.}
#'   \item{dup2}{Character, duplicated gene 2.}
#'   \item{type}{Factor of duplication types, with levels
#'               "SD" (segmental duplication),
#'               "TD" (tandem duplication), 
#'               "PD" (proximal duplication), and
#'               "DD" (dispersed duplication).}
#' }
#' @rdname get_tandem_proximal
#' @export
#' @examples
#' data(yeast_annot)
#' data(yeast_seq)
#' data(fungi_kaks)
#' scerevisiae_kaks <- fungi_kaks$saccharomyces_cerevisiae
#' 
#' # Get processed annotation for S. cerevisiae
#' pdata <- annotation <- syntenet::process_input(yeast_seq, yeast_annot)
#' annot <- pdata$annotation[[1]]
#' 
#' # Get duplicated pairs
#' pairs <- scerevisiae_kaks[, c("dup1", "dup2", "type")]
#' pairs$dup1 <- paste0("Sce_", pairs$dup1)
#' pairs$dup2 <- paste0("Sce_", pairs$dup2)
#' 
#' # Get tandem and proximal duplicates
#' td_pd_pairs <- get_tandem_proximal(pairs, annot)
#' 
get_tandem_proximal <- function(
        pairs = NULL, annotation_granges = NULL, proximal_max = 10
) {
    
    annot <- as.data.frame(annotation_granges)
    annot <- annot[, c("seqnames", "gene", "start", "end")]
    pairs$type <- as.character(pairs$type)
    ssd_pairs <- pairs[pairs$type == "DD", ]
    
    # Add chromosome number and order in the chromosome for each gene
    annot <- annot[order(annot$seqnames, annot$start), ]
    annot_bychr <- split(annot, annot$seqnames)
    annot_order <- Reduce(rbind, lapply(annot_bychr, function(x) {
        x$order <- seq_len(nrow(x))
        return(x[, c("seqnames", "gene", "order")])
    }))
    
    # Create df with `dup1`, `dup2`, `type`, `chr1`, `order1`, `chr2`, `order2`
    ssd_pos <- merge(ssd_pairs, annot_order, by.x = "dup1", by.y = "gene")
    names(ssd_pos)[c(4, 5)] <- c("chr_dup1", "order_dup1")
    ssd_pos <- merge(
        ssd_pos, annot_order, sort = FALSE, by.x = "dup2", by.y = "gene"
    )[, c(2, 1, 3, 4, 5, 6, 7)] # dup1, dup2, type, chr1, order1, chr2, order2
    names(ssd_pos)[c(6, 7)] <- c("chr_dup2", "order_dup2")
    
    # Find tandem and proximal duplicates
    ssd_pos$dist <- abs(ssd_pos$order_dup1 - ssd_pos$order_dup2)
    td_idx <- which(ssd_pos$chr_dup1 == ssd_pos$chr_dup2 & ssd_pos$dist == 1)
    pd_idx <- which(ssd_pos$chr_dup1 == ssd_pos$chr_dup2 & ssd_pos$dist > 1 &
                        ssd_pos$dist <= proximal_max)
    
    if(length(td_idx) > 0) { ssd_pos$type[td_idx] <- "TD" }
    if(length(pd_idx) > 0) { ssd_pos$type[pd_idx] <- "PD" }

    duplicates <- rbind(pairs[pairs$type != "DD", ], ssd_pos[, c(1, 2, 3)])
    l <- c("SD", "TD", "PD", "DD")
    duplicates$type <- factor(duplicates$type, levels = l)
    
    return(duplicates)
}


#' Get syntenic block ID for each gene in a gene pair
#' 
#' @param pair A 2-column data frame with gene pairs.
#' @param syn_df A 2-column data frame with gene ID in column 1, and synteny
#' block ID in column 2.
#'
#' @return A data frame of 4 columns as below:
#' \describe{
#'   \item{dup1}{Character, ID of duplicated gene 1.}
#'   \item{dup2}{Character, ID of duplicated gene 2.}
#'   \item{block1}{Numeric, syntenic block ID of gene 1.}
#'   \item{block2}{Numeric, syntenic block ID of gene 2.}
#' }
#'
#' @noRd
pairs_and_synblocks <- function(pairs, syn_df) {
    
    names(pairs)[1:2] <- c("dup1", "dup2")
    names(syn_df)[1:2] <- c("anchor", "block")
    
    pairs_ancestral <- merge(
        pairs[, c(1, 2)], syn_df, by.x = "dup1", by.y = "anchor",
        all.x = TRUE
    )
    pairs_ancestral <- merge(
        pairs_ancestral, syn_df, by.x = "dup2", by.y = "anchor",
        all.x = TRUE
    )
    names(pairs_ancestral)[c(3, 4)] <- c("block1", "block2")
    
    return(pairs_ancestral)
}


#' Classify gene pairs originating from transposon-derived duplications
#'
#' @param pairs A 3-column data frame with columns \strong{dup1}, \strong{dup2},
#' and \strong{type} indicating duplicated gene 1, duplicated gene 2, and
#' the mode of duplication associated with the pair. This data frame
#' is returned by \code{get_tandem_proximal()}.
#' @param blast_inter A list of data frames of length 1 
#' containing BLAST tabular output for the comparison between the target
#' species and an outgroup. Names of list elements must match the names of 
#' list elements in `annotation`. BLASTp, DIAMOND or simular programs must 
#' be run on processed sequence data as returned 
#' by \code{syntenet::process_input()}.
#' @param annotation A processed GRangesList or CompressedGRangesList object as
#' returned by \code{syntenet::process_input()}.
#' @param evalue Numeric scalar indicating the E-value threshold. 
#' Default: 1e-10.
#' @param anchors Numeric indicating the minimum required number of genes
#' to call a syntenic block, as in \code{syntenet::infer_syntenet}. 
#' Default: 5.
#' @param max_gaps Numeric indicating the number of upstream and downstream
#' genes to search for anchors, as in \code{syntenet::infer_syntenet}. 
#' Default: 25.
#' @param collinearity_dir Character indicating the path to the directory
#' where .collinearity files will be stored. If NULL, files will
#' be stored in a subdirectory of \code{tempdir()}. Default: NULL.
#' @param outgroup_coverage Numeric indicating the minimum percentage of 
#' outgroup species to use to consider genes as transposed duplicates. Only
#' valid if multiple outgroup species are present (see details below). Values
#' should range from 0 to 100. Default: 70.
#' 
#'
#' @return A 3-column data frame with the following variables:
#' \describe{
#'   \item{dup1}{Character, duplicated gene 1.}
#'   \item{dup2}{Character, duplicated gene 2.}
#'   \item{type}{Factor of duplication types, with levels
#'               "SD" (segmental duplication),
#'               "TD" (tandem duplication), 
#'               "PD" (proximal duplication), 
#'               "TRD" (transposon-derived duplication), and
#'               "DD" (dispersed duplication).}
#' }
#' 
#' @details 
#' If the list of interspecies DIAMOND tables contain comparisons of the
#' same species to multiple outgroups (e.g., 
#' 'speciesA_speciesB', 'speciesA_speciesC'), this function will check if
#' gene pairs are classified as transposed (i.e.,
#' only one gene is an ancestral locus) in each of the outgroup species,
#' and then calculate the percentage of outgroup species in which each pair
#' is considered 'transposed'. For instance, gene pair 1 is transposed based on
#' 30\% of the outgroup species, gene pair is considered as transposed based 
#' on  100\% of the outgroup species, gene pair 3 is considered as transposed
#' based on 0\% of the outgroup species, and so on. 
#' Parameter \strong{outgroup_coverage} lets you choose a minimum percentage 
#' cut-off to classify pairs as transposed.
#' 
#' @importFrom syntenet interspecies_synteny
#' @export
#' @rdname get_transposed
#' @examples 
#' data(diamond_inter)
#' data(diamond_intra)
#' data(yeast_seq)
#' data(yeast_annot)
#' data(fungi_kaks)
#' scerevisiae_kaks <- fungi_kaks$saccharomyces_cerevisiae
#' 
#' # Get processed annotation
#' pdata <- syntenet::process_input(yeast_seq, yeast_annot)
#' annotation <- pdata$annotation
#' 
#' # Get duplicated pairs
#' pairs <- scerevisiae_kaks[, c("dup1", "dup2", "type")]
#' pairs$dup1 <- paste0("Sce_", pairs$dup1)
#' pairs$dup2 <- paste0("Sce_", pairs$dup2)
#' 
#' # Classify pairs
#' trd <- get_transposed(pairs, diamond_inter, annotation)
#' 
#' annotation <- c(annotation, list(Cglabrata2 = annotation$Cglabrata))
#' blast_inter <- c(diamond_inter, list(Scerevisiae_Cglabrata2 = diamond_inter[[1]])) 
#' 
get_transposed <- function(
        pairs, blast_inter, annotation, 
        evalue = 1e-10, anchors = 5, max_gaps = 25,
        collinearity_dir = NULL, outgroup_coverage = 70
) {
    
    blast_inter <- lapply(blast_inter, function(x) return(x[x$evalue <= evalue, ]))
    pairs$type <- as.character(pairs$type)
    pairs_dd <- pairs[pairs$type == "DD", ]
    
    # Define directory where interspecies .collinearity files will be stored
    interdir <- collinearity_dir
    if(is.null(interdir)) {
        daytime <- format(Sys.time(), "%d_%b_%Y_%Hh%M")
        interdir <- file.path(tempdir(), paste0("inter_", daytime))
    }
    
    # Get name of target and outgroup species
    target <- unlist(lapply(names(annotation), function(x) {
        nfound <- sum(grepl(paste0(x, "_"), names(blast_inter)))
        found <- rep(x, nfound)
        return(found)
    }))
    
    outgroup <- unlist(lapply(names(annotation), function(x) {
        nfound <- sum(grepl(paste0(x, "$"), names(blast_inter)))
        found <- rep(x, nfound)
        return(found)
    }))

    # For each outgroup, get data frame indicating if pairs is tranposed
    trd_df <- lapply(seq_along(outgroup), function(n) {
        # Detect syntenic regions between `target` and `outgroup`
        syn <- syntenet::interspecies_synteny(
            blast_inter[n],
            annotation = annotation[c(target[n], outgroup[n])],
            inter_dir = interdir,
            anchors = anchors,
            max_gaps = max_gaps
        )
        
        # Read and parse interspecies synteny results
        parsed_syn <- collinearity2blocks(syn)[, c("anchor2", "block")]
        parsed_syn <- parsed_syn[!duplicated(parsed_syn$anchor2), ]
        
        pairs_ancestral <- pairs
        pairs_ancestral$ancestral <- FALSE
        if(!is.null(parsed_syn)) {
            
            # Find TRD-derived genes (only one member of pair in ancestral loci)
            pairs_ancestral <- pairs_and_synblocks(pairs_dd, parsed_syn)
            
            nas <- apply(pairs_ancestral[, c(3, 4)], 1, function(x) return(sum(is.na(x))))
            pairs_ancestral$ancestral <- ifelse(nas == 1, TRUE, FALSE)
            pairs_ancestral <- pairs_ancestral[, c("dup1", "dup2", "ancestral")]
        }
        
        return(pairs_ancestral)
    })
    nout <- length(trd_df)
    trd_df <- Reduce(function(x, y) merge(x, y, by = c("dup1", "dup2"), all = TRUE), trd_df)
    names(trd_df)[seq(3, nout+2, 1)] <- paste0("ancestral", seq_len(nout))
    
    # Calculate percentage of outgroups in which pair is classified as transposed
    perc_trd <- (rowSums(trd_df[, -c(1,2), drop = FALSE]) / nout) * 100
    
    # Classify pairs as TRD if percentage >= outgroup_coverage
    trd_df$type <- ifelse(perc_trd >= outgroup_coverage, "TRD", "DD")
    final <- rbind(
        pairs[pairs$type != "DD", ],
        trd_df[, c("dup1", "dup2", "type")]
    )
    final$type <- factor(final$type, levels = c("SD", "TD", "PD", "TRD", "DD"))
    
    return(final)
}



#' Classify TRD genes as derived from either DNA transposons or retrotransposons
#'
#' @param pairs A 3-column data frame with columns \strong{dup1}, \strong{dup2},
#' and \strong{type} indicating duplicated gene 1, duplicated gene 2, and
#' the mode of duplication associated with the pair. This data frame
#' is returned by \code{get_transposed()}.
#' @param intron_counts A 2-column data frame with columns \strong{gene}
#' and \strong{introns} indicating the number of introns for each gene,
#' as returned by \code{get_intron_counts}.
#'
#' @return A 3-column data frame with the following variables:
#' \describe{
#'   \item{dup1}{Character, duplicated gene 1.}
#'   \item{dup2}{Character, duplicated gene 2.}
#'   \item{type}{Factor of duplication types, with levels
#'               "SD" (segmental duplication),
#'               "TD" (tandem duplication), 
#'               "PD" (proximal duplication), 
#'               "dTRD" (DNA transposon-derived duplication),
#'               "rTRD" (retrotransposon-derived duplication), and
#'               "DD" (dispersed duplication).}
#' }
#' 
#' @rdname get_transposed_classes
#' @export
#' @examples
#' data(diamond_inter)
#' data(diamond_intra)
#' data(yeast_seq)
#' data(yeast_annot)
#' data(fungi_kaks)
#' scerevisiae_kaks <- fungi_kaks$saccharomyces_cerevisiae
#' 
#' # Get processed annotation
#' pdata <- syntenet::process_input(yeast_seq, yeast_annot)
#' annotation <- pdata$annotation
#' 
#' # Get duplicated pairs
#' pairs <- scerevisiae_kaks[, c("dup1", "dup2", "type")]
#' pairs$dup1 <- paste0("Sce_", pairs$dup1)
#' pairs$dup2 <- paste0("Sce_", pairs$dup2)
#' 
#' # Classify pairs
#' trd <- get_transposed(pairs, diamond_inter, annotation)
#' 
#' # Create TxDb object from GRanges
#' library(txdbmaker)
#' txdb <- txdbmaker::makeTxDbFromGRanges(yeast_annot[[1]])
#'
#' # Get intron counts
#' intron_counts <- get_intron_counts(txdb)
#'
#' # Get TRD classes
#' trd_classes <- get_transposed_classes(trd, intron_counts)
#'
get_transposed_classes <- function(pairs, intron_counts) {
    
    # Get TRD pairs
    pairs$type <- as.character(pairs$type)
    tpairs <- pairs[pairs$type == "TRD", ]
    
    final_pairs <- pairs
    if(nrow(tpairs) > 0) {
        
        id <- unique(gsub("_.*", "", tpairs$dup1))
        tpairs$dup1 <- gsub("^[a-zA-Z]{2,5}_", "", tpairs$dup1)
        tpairs$dup2 <- gsub("^[a-zA-Z]{2,5}_", "", tpairs$dup2)
        
        # Combine `tpairs` and `intron_counts`
        pairs_ic <- merge(
            tpairs, intron_counts, by.x = "dup1", by.y = "gene", all.x = TRUE
        )
        pairs_ic <- merge(
            pairs_ic, intron_counts, by.x = "dup2", by.y = "gene", all.x = TRUE
        )
        names(pairs_ic)[c(4,5)] <- c("introns1", "introns2")
        
        # Create a column with number of genes in pair with 0 introns
        pairs_ic$count <- apply(pairs_ic[, 4:5], 1, function(x) sum(x == 0))
        
        # 'rTRD' if only one gene has no introns, as 'dTRD' otherwise
        pairs_ic$type <- ifelse(pairs_ic$count == 1, "rTRD", "dTRD")
        
        final_pairs <- pairs_ic[, c("dup1", "dup2", "type")]
        
        # Add species IDs back
        final_pairs$dup1 <- paste0(id, "_", final_pairs$dup1)
        final_pairs$dup2 <- paste0(id, "_", final_pairs$dup2)
        
        final_pairs <- rbind(pairs[pairs$type != "TRD", ], final_pairs)
    }
    
    l <- c("SD", "TD", "PD", "rTRD", "dTRD", "DD")
    final_pairs$type <- factor(final_pairs$type, levels = l)
    
    return(final_pairs)
}

