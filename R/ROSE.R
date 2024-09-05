#' Fetches TSSes or promoters from a TxDb object
#'
#' @param TxDb TxDb object to fetch TSSes from.
#' @param expansion.distance Integer __distance__ to add up/downstream of TSSes.
#' @return GRanges object of TSSes for a TxDb object.
#'
#' @export
#' @importFrom GenomicRanges promoters
#' @author Jacqueline Myers
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' tsses <- fetch_tss(TxDb, expansion.distance = 1000)
fetch_tss <- function(TxDb, expansion.distance = 0) {

    if (expansion.distance > 0) {
        y <- promoters(TxDb, upstream = expansion.distance, downstream =  expansion.distance)
    } else {
        y <- promoters(TxDb)
    }
    return(y)
}


#' Normalizes the coverage of a candidate enhancer based on total reads and region size
#'
#' This replicates the normalization scheme from ROSE directly (sorta, in theory).
#'
#' @details
#' First it will extract the total number of reads in the bam file to calculate a scaling factor (total mapped reads/1M).
#' It will calculate region length, required for region density = (total reads in region/ region length).
#' Finally, it will compute normalized signal = read density / scaling factor.
#'
#'
#' @param candidate.enhancer RangedSummarizedExperiment object containing counts for ChIP and input.
#' @param stitched.regions GRanges object containing stitched regions.
#' @return RangedSummarizedExperiment object with normalized read density as a new assay - NormReadDensity.
#'
#' @importFrom GenomicRanges width
#' @importFrom BiocParallel bplapply
#' @export
#'
#' @author Jacqueline Myers
enhancer_norm <- function(candidate.enhancer, stitched.regions) {
    # compute the scaling factor
    # sf returns a sf for the ChIP enriched and the input
    sf <- candidate.enhancer$totals/1000000
    # calculate region sizes and save it in obj for later use
    candidate.enhancer@assays@data$counts <- cbind(candidate.enhancer@assays@data$counts,width(stitched.regions))
    # loop through each region and compute normalized read density
    # use lapply and bplapply
    candidate.enhancer@assays@data$NormReadDensity <- apply(candidate.enhancer@assays@data$counts,1, function(x) {
        if (length(x) == 3) {
            # x[1] = chip counts
            # X[2] = input counts
            # x[3] = region size
            # sf[1] = scaling factor for chip
            # sf[2] = scaling factor for input
            ((x[1] / x[3]) / sf[1]) - ((x[2] / x[3]) / sf[2])
        } else {
            # x[1] = chip counts
            # X[2] = region size
            # sf[1] = scaling factor for chip
            ((x[1] / x[2]) / sf[1])
        }
    })
    return(candidate.enhancer)
}
