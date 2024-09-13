#' Fetches TSSes or promoters from a TxDb object
#'
#' This function fetches TSS regions and optionally expands them by a specified number of bases.
#'
#' @param TxDb TxDb object to fetch TSSes from.
#' @param expansion.distance Number of bases to add up/downstream of TSSes.
#' @return GRanges object of TSSes or promoters for each transcrip in `TxDb`.
#'
#' @export
#' @importFrom GenomicRanges promoters
#' 
#' @author Jacqueline Myers, Jared Andrews
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' tsses <- fetch_tss(TxDb, expansion.distance = 1000)
fetch_tss <- function(TxDb, expansion.distance = 0) {
    if (expansion.distance > 0) {
        y <- promoters(TxDb, upstream = expansion.distance, downstream = expansion.distance)
    } else {
        y <- promoters(TxDb, upstream = 0, downstream = 0)
    }
    y
}


#' Extend reads directionally upstream and downstream
#'
#' Extends reads by a specified number of bases upstream or downstream
#' based on strand orientation.
#'
#' @param regions A GRanges object containing genomic regions.
#' @param upstream Number of bases to extend upstream. Default is 0.
#' @param downstream Number of bases to extend downstream. Default is 0.
#'
#' @return A GRanges object with extended regions.
#'
#' @export
#' @importFrom GenomicRanges strand start end ranges
#' @importFrom IRanges IRanges
#' 
#' @author Jared Andrews
#'
#' @examples
#' regions <- GRanges(seqnames = "chr1", ranges = IRanges(100, 200), strand = "+")
#' extended <- extend_reads(regions, upstream = 50, downstream = 100)
extend_reads <- function(regions, upstream = 0, downstream = 0) {
    pos <- strand(regions) == "+" | strand(regions) == "*"
    ext_start <- start(regions) - ifelse(pos, upstream, downstream)
    ext_end <- end(regions) + ifelse(pos, downstream, upstream)
    ranges(regions) <- IRanges(ext_start, ext_end)
    regions
}
