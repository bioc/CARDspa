#' Extension of CARD into a reference-free version of
#' deconvolution: CARDfree.
#'
#' @param markerlist a list of marker genes, with each element of the list
#' being the vector of cell type specific marker genes
#' @param spatial_count Raw spatial resolved transcriptomics data, each column
#' is a spatial location, and each row is a gene.
#' @param spatial_location data frame, with two columns representing the x and
#' y coordinates of the spatial location. The rownames of this data frame
#' should match eaxctly with the columns of the spatial_count.
#' @param spe a \code{SpatialExperiment} object containing spatial data in 
#' the \code{counts} assay, and spatial coordinates in the spatialCoords.
#' @param mincountgene Minimum counts for each gene
#' @param mincountspot Minimum counts for each spatial location
#' @importFrom Rcpp sourceCpp
#' @importFrom fields rdist
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' @return Returns a SpatialExperiment object with estimated cell type 
#' proportion stored in object$Proportion_CARD. Because this is a 
#' reference-free version, the columns of estimated proportion is not cell 
#' type but cell type cluster
#'
#' @export
#' @examples
#' library(NMF)
#' library(RcppArmadillo)
#' data(markerList)
#' data(spatial_count)
#' data(spatial_location)
#' CARDfree_obj <- CARD_refFree(
#' markerlist = markerList[8:16],
#' spatial_count = spatial_count[1:2500, ],
#' spatial_location = spatial_location,
#' mincountgene = 100,
#' mincountspot = 5
#' )
#'
CARD_refFree <- function(
        markerlist, 
        spatial_count, 
        spatial_location,
        mincountgene = 100, 
        mincountspot = 5,
        spe = NULL) {
    ##### Check
    if(is.null(spatial_count) || is.null(spatial_location)){
        if(is.null(spe)){
            stop("Please provide SpatialExperiment object")}
    }else{
        if (!is.vector(spatial_count) && 
            !is.matrix(spatial_count) && 
            !inherits(spatial_count, "sparseMatrix")) {
            stop("SRT counts has to be of following forms: vector,
        matrix or sparseMatrix")
        }
    }
    if(is.null(markerlist)){
        stop("Please provide the marker list for CARDfree")
    }
    
    ##### Create CARD refFree obejct
    CARDfree_object <- createCARDfreeObject(
        markerlist = markerlist,
        spatial_count = spatial_count,
        spatial_location = spatial_location,
        mincountgene = mincountgene,
        mincountspot = mincountspot,
        spe = spe
    )
    
    #### Run CARD refFree deconvolution
    ##### load in spatial transcriptomics data stored in CARDfree_object
    spatial_countmat <- CARDfree_object@spatial_countMat
    spatial_location <- CARDfree_object@spatial_location
    
    ##### load in markerlist number of cell type clusters
    markerlist <- CARDfree_object@markerList
    numK <- length(markerlist)
    marker <- unique(unlist(markerlist))
    message(
        "## Number of unique marker genes: ", 
        length(marker), 
        " for ",
        numK, 
        " cell types ...\n"
    )
    commongene <- intersect(
        toupper(rownames(spatial_countmat)),
        toupper(marker)
    )
    
    ##### remove mitochondrial and ribosomal genes
    commongene <- commongene[!(commongene %in% 
                            commongene[grep("mt-",commongene)])]
    if (length(commongene) < numK * 10) {
        stop("## STOP! The average number of unique marker genes for each cell
            type is less than 20 ...\n")
    }
    xinput <- spatial_countmat[order(rownames(spatial_countmat)), ]
    xinput <- xinput[order(rownames(xinput)), ]
    xinput <- xinput[toupper(rownames(xinput)) %in% commongene, ]
    xinput <- xinput[rowSums(xinput) > 0, ]
    xinput <- xinput[, colSums(xinput) > 0]
    xinput_norm <- sweep(xinput, 2, colSums(xinput), "/")
    
    ##### initialization
    if (ncol(xinput_norm) < 5000) {
        if (.Platform$OS.type == "windows") {
            sink("NUL")
        } else {
            sink("/dev/null")
        }
        nmfout <- invisible(NMF::nmf(as.matrix(xinput_norm), numK))
        sink()
        B <- nmfout@fit@W
        vint1 <- as.matrix(t(nmfout@fit@H))
        rownames(vint1) <- colnames(xinput_norm)
    } else {
        if (.Platform$OS.type == "windows") {
            sink("NUL")
        } else {
            sink("/dev/null")
        }
        if (requireNamespace("RcppML", quietly = TRUE)) {
            nmfout <- invisible(RcppML::nmf(as.matrix(xinput_norm), numK))
            sink()
            B <- nmfout$w
            rownames(B) <- rownames(xinput_norm)
            vint1 <- as.matrix(t(nmfout$h))
            rownames(vint1) <- colnames(xinput_norm)
        } else {
            # fallback if RcppML not available
            nmfout <- invisible(NMF::nmf(as.matrix(xinput_norm), rank = numK))
            sink()
            B <- nmfout@fit@W
            vint1 <- as.matrix(t(nmfout@fit@H))
            rownames(vint1) <- colnames(xinput_norm)
        }
    }
    spatial_location <- spatial_location[rownames(spatial_location) %in%
                                        colnames(xinput_norm), ]
    spatial_location <- spatial_location[
        match(colnames(xinput_norm),rownames(spatial_location)), ]
    norm_cords <- spatial_location[, c("x", "y")]
    norm_cords$x <- norm_cords$x - min(norm_cords$x)
    norm_cords$y <- norm_cords$y - min(norm_cords$y)
    scalefactor <- max(norm_cords$x, norm_cords$y)
    norm_cords$x <- norm_cords$x / scalefactor
    norm_cords$y <- norm_cords$y / scalefactor
    
    ##### Euclidiean distance
    ed <- rdist(as.matrix(norm_cords)) ## Euclidean distance matrix
    b <- rep(0, ncol(B))
    
    ##### parameters that need to be set
    isigma <- 0.1
    epsilon <- 1e-04 #### convergence epsion
    phi <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99) #### grided values for phi
    kernel_mat <- exp(-ed^2 / (2 * isigma^2))
    diag(kernel_mat) <- 0
    rm(ed)
    rm(xinput)
    rm(norm_cords)
    gc()
    
    ##### scale the xinput_norm and B to speed up the convergence.
    mean_X <- mean(xinput_norm)
    mean_B <- mean(B)
    xinput_norm <- xinput_norm * 0.1 / mean_X
    B <- B * 0.1 / mean_B
    gc()
    reslist <- list()
    obj <- c()
    for (iphi in seq_len(length(phi))) {
        res <- CARDfree(
            XinputIn = as.matrix(xinput_norm),
            UIn = as.matrix(B),
            WIn = kernel_mat,
            phiIn = phi[iphi],
            max_iterIn = 1000,
            epsilonIn = epsilon,
            initV = vint1,
            initb = rep(0, ncol(B)),
            initSigma_e2 = 0.1,
            initLambda = rep(10, ncol(B))
        )
        rownames(res$V) <- colnames(xinput_norm)
        colnames(res$V) <- paste0("CT", seq_len(ncol(B)))
        reslist[[iphi]] <- res
        obj <- c(obj, res$Obj)
    }
    optimal_ind <- which(obj == max(obj))
    optimal <- optimal_ind[length(optimal_ind)]
    optimalphi <- phi[optimal]
    optimalres <- reslist[[optimal]]
    message("## Deconvolution Finish! ...\n")
    CARDfree_object@info_parameters$phi <- optimalphi
    CARDfree_object@Proportion_CARD <- sweep(
        optimalres$V, 
        1,
        rowSums(optimalres$V),
        "/"
    )
    CARDfree_object@algorithm_matrix <- list(
        B = optimalres$B * mean_B / 0.1,
        Xinput_norm = 
            xinput_norm * mean_X / 0.1, 
        Res = optimalres
    )
    CARDfree_object@spatial_location <- spatial_location
    CARDfree_object@estimated_refMatrix <- optimalres$B * mean_B / 0.1
    rownames(CARDfree_object@estimated_refMatrix) <- rownames(B)
    colnames(CARDfree_object@estimated_refMatrix) <- colnames(B)
    
    ##### store as a SpatialExperiemnt object
    spe_out <- SpatialExperiment(
        assay = list(spatial_countMat = CARDfree_object@spatial_countMat), 
        spatialCoords = as.matrix(CARDfree_object@spatial_location))
    ## remove sample_id 
    spe_out@colData@listData[["sample_id"]] <- NULL
    colData(spe_out)$Proportion_CARD <- CARDfree_object@Proportion_CARD
    metadata(spe_out)$project <- CARDfree_object@project
    metadata(spe_out)$info_parameters <- CARDfree_object@info_parameters
    metadata(spe_out)$algorithm_matrix <- CARDfree_object@algorithm_matrix
    metadata(spe_out)$markerList <- CARDfree_object@markerList
    metadata(spe_out)$estimated_refMatrix <- 
        CARDfree_object@estimated_refMatrix
    return(spe_out)
}
