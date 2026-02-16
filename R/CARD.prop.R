#' Construct the mean gene expression basis matrix (B), this is the
#' faster version
#'
#' @param sc_eset S4 class for storing data from single-cell experiments. 
#' This format is usually created by the package SingleCellExperiment with 
#' stored counts, along with the usual metadata for genes and cells.
#' @param ct_select vector of cell type names that you are interested in to
#' deconvolute, default as NULL. If NULL, then use all cell types provided by
#' single cell dataset;
#' @param ct_varname character, the name of the column in metaData that
#' specifies the cell type annotation information
#' @param sample_varname character,the name of the column in metaData that
#' specifies the sample information. If NULL, we just use the whole as one
#' sample.
#'
#' @importFrom SummarizedExperiment assays colData
#' @importFrom wrMisc rowGrpMeans
#' @importMethodsFrom SummarizedExperiment colData<-
#' @importFrom stats aggregate reshape
#' @return Return a list of basis (B) matrix
#'

create_ref <- function(
        sc_eset, 
        ct_select = NULL, 
        ct_varname,
        sample_varname = NULL) {
    ##### Prepare
    if (is.null(ct_select)) {
        ct_select <- unique(colData(sc_eset)[, ct_varname])
    }
    ct_select <- ct_select[!is.na(ct_select)]
    count_mat <- as(assays(sc_eset)$counts, "sparseMatrix")
    ct_id <- droplevels(as.factor(colData(sc_eset)[, ct_varname]))
    if (is.null(sample_varname)) {
        colData(sc_eset)$sampleID <- "Sample"
        sample_varname <- "sampleID"
    }
    sample_id <- as.character(colData(sc_eset)[, sample_varname])
    ct_sample_id <- paste(ct_id, sample_id, sep = "$*$")
    colsums_count_mat <- colSums(count_mat)
    colsums_count_mat_ct <- aggregate(
        colsums_count_mat ~ ct_id + sample_id,
        FUN = "sum"
    )
    colsums_count_mat_ct_wide <- reshape(
        colsums_count_mat_ct,
        idvar = "sample_id",
        timevar = "ct_id",
        direction = "wide"
    )
    colnames(colsums_count_mat_ct_wide) <-
        gsub("colsums_count_mat.", "", colnames(colsums_count_mat_ct_wide))
    rownames(colsums_count_mat_ct_wide) <- colsums_count_mat_ct_wide$sample_id
    colsums_count_mat_ct_wide$sample_id <- NULL
    tbl <- table(sample_id, ct_id)
    colsums_count_mat_ct_wide <-
        colsums_count_mat_ct_wide[,
            match(colnames(tbl),colnames(colsums_count_mat_ct_wide))]
    colsums_count_mat_ct_wide <-
        colsums_count_mat_ct_wide[
            match(rownames(tbl),rownames(colsums_count_mat_ct_wide)), ]
    
    ##### Calculate gene expression basis
    s_jk <- colsums_count_mat_ct_wide / tbl
    s_jk <- as.matrix(s_jk)
    s_jk[s_jk == 0] <- NA
    s_jk[!is.finite(s_jk)] <- NA
    S <- colMeans(s_jk, na.rm = TRUE)
    S <- S[match(unique(ct_id), names(S))]
    if (nrow(count_mat) > 10000 & ncol(count_mat) > 50000) {
        ### to save memory
        seqID <- seq(1, nrow(count_mat), by = 10000)
        theta_s_rowmean <- NULL
        for (igs in seqID) {
            if (igs != seqID[length(seqID)]) {
                theta_s_rowmean_tmp <- wrMisc::rowGrpMeans(
                    as.matrix(count_mat[(igs:(igs + 9999)), ]),
                    grp = ct_sample_id,
                    na.rm = TRUE
                )
            } else {
                theta_s_rowmean_tmp <- wrMisc::rowGrpMeans(
                    as.matrix(count_mat[igs:nrow(count_mat), ]),
                    grp = ct_sample_id,
                    na.rm = TRUE
                )
            }
            theta_s_rowmean <- rbind(theta_s_rowmean, theta_s_rowmean_tmp)
        }
    } else {
        theta_s_rowmean <- wrMisc::rowGrpMeans(as.matrix(count_mat),
            grp = ct_sample_id,
            na.rm = TRUE
        )
    }
    tbl_sample <- table(ct_sample_id)
    tbl_sample <- tbl_sample[match(
        colnames(theta_s_rowmean),
        names(tbl_sample))]
    theta_s_rowsums <- sweep(theta_s_rowmean, 2, tbl_sample, "*")
    theta_s <- sweep(theta_s_rowsums, 2, colSums(theta_s_rowsums), "/")
    grp <- vapply(
        strsplit(colnames(theta_s), split = "$*$", fixed = TRUE),
        `[`,
        FUN.VALUE = character(1),
        1
    )
    theta <- wrMisc::rowGrpMeans(theta_s, grp = grp, na.rm = TRUE)
    theta <- theta[, match(unique(ct_id), colnames(theta))]
    S <- S[match(colnames(theta), names(S))]
    basis <- sweep(theta, 2, S, "*")
    colnames(basis) <- colnames(theta)
    rownames(basis) <- rownames(theta)
    return(list(basis = basis))
}

#' Select Informative Genes used in the deconvolution
#'
#' @param basis Reference basis matrix.
#' @param sc_eset scRNAseq data along with meta data stored in the S4 class
#' format (SingleCellExperiment).
#' @param commongene common genes between scRNAseq count data and spatial
#' resolved transcriptomics data.
#' @param ct_select vector of cell type names that you are interested in to
#' deconvolute, default as NULL. If NULL, then use all cell types provided by
#' single cell dataset;
#' @param ct_varname character, the name of the column in metaData that
#' specifies the cell type annotation information
#'
#' @importFrom SummarizedExperiment assays colData
#' @importFrom stats quantile var
#' @return a vector of informative genes selected
#'

select_info <- function(basis, sc_eset, commongene, ct_select, ct_varname) {
    ##### log2 mean fold change >0.5
    gene1 <- lapply(ct_select, function(ict) {
        rest <- rowMeans(basis[, colnames(basis) != ict])
        fc <- log((basis[, ict] + 1e-06)) - log((rest + 1e-06))
        rownames(basis)[fc > 1.25 & basis[, ict] > 0]
    })
    gene1 <- unique(unlist(gene1))
    gene_int <- intersect(gene1, commongene)
    counts <- assays(sc_eset)$counts
    counts <- counts[rownames(counts) %in% gene_int, ]
    
    ##### only check the cell type that contains at least 2 cells
    ct_select <- names(table(colData(sc_eset)[, ct_varname])
        )[table(colData(sc_eset)[, ct_varname]) > 1]
    sd_within <- vapply(ct_select, function(ict) {
        temp <- counts[, colData(sc_eset)[, ct_varname] == ict]
        apply(temp, 1, stats::var) / apply(temp, 1, mean)
    }, FUN.VALUE = numeric(nrow(counts)))
    
    ##### remove the outliers that have high dispersion across cell types
    gene2 <- rownames(sd_within)[apply(sd_within,1,mean,na.rm = TRUE) < 
                                quantile(
                                    apply(sd_within,1,mean,na.rm = TRUE),
                                    prob = 0.99,
                                    na.rm = TRUE)]
    return(gene2)
}


#' Spatially Informed Cell Type Deconvolution for Spatial Transcriptomics
#' by CARD
#'
#' @param sc_count Raw scRNA-seq count data, each column is a cell and each row
#' is a gene.
#' @param sc_meta data frame, with each row representing the cell type and/or
#' sample information of a specific cell. The row names of this data frame
#' should match exactly with the column names of the sc_count data
#' @param spatial_count Raw spatial resolved transcriptomics data, each column
#' is a spatial location, and each row is a gene.
#' @param spatial_location data frame, with two columns representing the x and
#' y coordinates of the spatial location. The rownames of this data frame
#' should match eaxctly with the columns of the spatial_count.
#' @param sce a \code{SingleCellExperiment} object containing scRNA-seq count 
#' data in the \code{counts} assay, and cell types and sample information in 
#' the colData.
#' @param spe a \code{SpatialExperiment} object containing spatial data in 
#' the \code{counts} assay, and spatial coordinates in the spatialCoords.
#' @param ct_varname character, the name of the column in metaData that
#' specifies the cell type annotation information
#' @param ct_select vector of cell type names that you are interested in to
#' deconvolute, default as NULL. If NULL, then use all cell types provided by
#' single cell dataset;
#' @param sample_varname character,the name of the column in metaData that
#' specifies the sample information. If NULL, we just use the whole as one
#' sample.
#' @param mincountgene Minimum counts for each gene
#' @param mincountspot Minimum counts for each spatial location
#' @importFrom Rcpp sourceCpp
#' @importFrom MCMCpack rdirichlet
#' @importFrom fields rdist
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' @return Returns a SpatialExperiment object with estimated cell type 
#' proportion stored in object$Proportion_CARD.
#'
#' @export
#' @examples
#' library(NMF)
#' library(RcppArmadillo)
#' data(spatial_count)
#' data(spatial_location)
#' data(sc_count)
#' data(sc_meta)
#' CARD_obj <- CARD_deconvolution(
#'     sc_count = sc_count,
#'     sc_meta = sc_meta,
#'     spatial_count = spatial_count,
#'     spatial_location = spatial_location,
#'     ct_varname = "cellType",
#'     ct_select = unique(sc_meta$cellType),
#'     sample_varname = "sampleInfo",
#'     mincountgene = 100,
#'     mincountspot = 5
#' )
CARD_deconvolution <- function(
        sc_count, 
        sc_meta, 
        spatial_count, 
        spatial_location,
        ct_varname, 
        ct_select, 
        sample_varname,
        mincountgene = 100, 
        mincountspot = 5,
        sce = NULL, 
        spe = NULL) {
    ##### Check
    if(is.null(sc_count) ||
        is.null(sc_meta) || 
        is.null(spatial_count) || 
        is.null(spatial_location)){
        if(is.null(sce) || is.null(spe)){
            stop("Please provide SingleCellExperiment object and 
            SpatialExperiment object")
        }
    }else{
        if (!is.vector(sc_count) && 
            !is.matrix(sc_count) && 
            !inherits(sc_count, "sparseMatrix")) {
            stop("scRNASeq counts has to be of following forms: vector,
        matrix or sparseMatrix")
        }
        if (!is.vector(spatial_count) && 
            !is.matrix(spatial_count) && 
            !inherits(spatial_count, "sparseMatrix")) {
            stop("SRT counts has to be of following forms: vector,
        matrix or sparseMatrix")
        }
    }
    if (is.null(ct_varname)) {
        stop("Please provide the column name indicating the cell type 
        information in the meta data of scRNA-seq")
    }
    if(!is.null(sample_varname) && 
        !is.character(sample_varname) || 
        length(sample_varname) > 1){
        stop("Please input valid sample variable name")
    }
    if(!is.null(ct_select) && 
        !is.character(ct_select) ){
        stop("Please input valid cell types")
    }
    
    ##### create CARD object
    CARD_object <- createCARDObject(
        sc_count = sc_count,
        sc_meta = sc_meta,
        sce = sce,
        spe = spe,
        spatial_count = spatial_count,
        spatial_location = spatial_location,
        ct_varname = ct_varname,
        ct_select = ct_select,
        sample_varname = sample_varname,
        mincountgene = mincountgene,
        mincountspot = mincountspot) 
    
    ##### run deconvolution
    ct_select <- CARD_object@info_parameters$ct.select
    ct_varname <- CARD_object@info_parameters$ct.varname
    sample_varname <- CARD_object@info_parameters$sample.varname
    message("## create reference matrix from scRNASeq...\n")
    sc_eset <- CARD_object@sc_eset
    basis_ref <- create_ref(sc_eset, ct_select, ct_varname, sample_varname)
    basis <- basis_ref$basis
    basis <- basis[, colnames(basis) %in% ct_select]
    basis <- basis[, match(ct_select, colnames(basis))]
    spatial_count <- CARD_object@spatial_countMat
    commongene <- intersect(rownames(spatial_count), rownames(basis))
    
    ##### remove mitochondrial and ribosomal genes
    commongene <- commongene[!(commongene %in%
                    commongene[grep("mt-", commongene)])]
    message("## Select Informative Genes! ...\n")
    common <- select_info(basis, sc_eset, commongene, ct_select, ct_varname)
    xinput <- spatial_count
    rm(spatial_count)
    B <- basis
    rm(basis)
    
    ##### match the common gene names
    xinput <- xinput[order(rownames(xinput)), ]
    B <- B[order(rownames(B)), ]
    B <- B[rownames(B) %in% common, ]
    xinput <- xinput[rownames(xinput) %in% common, ]
    
    ##### filter out non expressed genes or cells again
    xinput <- xinput[rowSums(xinput) > 0, ]
    xinput <- xinput[, colSums(xinput) > 0]
    
    ##### normalize count data
    colsumvec <- colSums(xinput)
    xinput_norm <- sweep(xinput, 2, colsumvec, "/")
    B <- B[rownames(B) %in% rownames(xinput_norm), ]
    B <- B[match(rownames(xinput_norm), rownames(B)), ]
    
    ##### spatial location
    spatial_location <- CARD_object@spatial_location
    spatial_location <- spatial_location[rownames(spatial_location) %in%
                                    colnames(xinput_norm), ]
    spatial_location <- spatial_location[match(
        colnames(xinput_norm),
        rownames(spatial_location)
    ), ]
    
    #### normalize the coordinates without changing the shape and relative
    ##### position
    norm_cords <- spatial_location[, c("x", "y")]
    norm_cords$x <- norm_cords$x - min(norm_cords$x)
    norm_cords$y <- norm_cords$y - min(norm_cords$y)
    scalefactor <- max(norm_cords$x, norm_cords$y)
    norm_cords$x <- norm_cords$x / scalefactor
    norm_cords$y <- norm_cords$y / scalefactor
    
    ##### initialize the proportion matrix
    ed <- rdist(as.matrix(norm_cords)) ## Euclidean distance matrix
    message("## Deconvolution Starts! ...\n")
    vint1 <- as.matrix(rdirichlet(ncol(xinput_norm), rep(10, ncol(B))))
    colnames(vint1) <- colnames(B)
    rownames(vint1) <- colnames(xinput_norm)
    b <- rep(0, length(ct_select))
    
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
        res <- CARDref(
            XinputIn = as.matrix(xinput_norm), 
            UIn = as.matrix(B),
            WIn = kernel_mat,
            phiIn = phi[iphi],
            max_iterIn = 1000,
            epsilonIn = epsilon,
            initV = vint1,
            initb = rep(0, ncol(B)),
            initSigma_e2 = 0.1,
            initLambda = rep(10, length(ct_select))
        )
        rownames(res$V) <- colnames(xinput_norm)
        colnames(res$V) <- colnames(B)
        reslist[[iphi]] <- res
        obj <- c(obj, res$Obj)
    }
    optimal_ind <- which(obj == max(obj))
    optimal <- optimal_ind[length(optimal_ind)]
    optimalphi <- phi[optimal]
    optimalres <- reslist[[optimal]]
    message("## Deconvolution Finish! ...\n")
    CARD_object@info_parameters$phi <- optimalphi
    CARD_object@Proportion_CARD <- sweep(
        optimalres$V, 1, rowSums(optimalres$V),
        "/"
    )
    CARD_object@algorithm_matrix <- list(
        B = B * mean_B / 0.1,
        Xinput_norm = xinput_norm * mean_X / 0.1, 
        Res = optimalres
    )
    CARD_object@spatial_location <- spatial_location
    
    ##### return SpatialExperiment Object
    spe_out <- SpatialExperiment(
        assay = list(spatial_countMat = CARD_object@spatial_countMat), 
        spatialCoords = as.matrix(CARD_object@spatial_location)
    )
    colData(spe_out)$Proportion_CARD <- CARD_object@Proportion_CARD
    
    ##### remove sample_id 
    spe_out@colData@listData[["sample_id"]] <- NULL
    metadata(spe_out)$project <- CARD_object@project
    metadata(spe_out)$info_parameters <- CARD_object@info_parameters
    metadata(spe_out)$algorithm_matrix <- CARD_object@algorithm_matrix
    metadata(spe_out)$sc_eset <- CARD_object@sc_eset
    return(spe_out)
}
