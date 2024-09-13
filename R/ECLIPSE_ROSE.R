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


#' Get region signal from BAM files
#'
#' Calculates the libary size-normalized coverage of reads within specified regions from a 
#' signal BAM file with optional subtraction of a control BAM file.
#' Allows read extension and additional processing like setting negative coverage to zero.
#' 
#' @details The defaults of this function are set to so as to match the functionality
#' of ROSE as closely as possible.
#' 
#' To detail the process:
#'  - The total number of reads in the signal BAM file is calculated.
#'  - The signal reads are extended downstream by a specified number of bases (200 bp by default).
#'  - The coverage for each basepair in the regions of interest are calculated. 
#'      ROSE does this manually by calling samtools for each region, which is slow.
#'  - Basepairs with coverage below a specified threshold (`floor`, 1 by default) are removed.
#'  - For each region, the coverage is summed and divided by the total number of reads to get the signal.
#'
#' @param sample.bam Character or BamFile object representing the signal BAM file.
#' @param regions Character or GRanges object representing genomic regions of interest.
#' @param control.bam Character or BamFile object representing the control BAM file.
#'   Default is `NULL`.
#' @param floor Numeric value for the minimum coverage threshold to count for region.
#'   Default is 0.
#' @param read.ext Numeric value for extending reads downstream. Default is 200.
#'
#' @return A GRanges object for `regions` with additional columns for sample and control (if provided) signal.
#'   The metadata of the object will also contain the scaling factor for library size normalization for the `sample.bam` 
#'   and `control.bam` (if provided).
#'
#' @export
#' @importFrom Rsamtools BamFile idxstatsBam
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges granges
#' @importFrom HelloRanges do_bedtools_coverage
#' @importFrom IRanges IRanges
#' @importFrom genomation readBed
#'
#' @author Jared Andrews
#'
#' @examples
#' sample.bam <- "path/to/sample.bam"
#' regions <- GRanges(seqnames = "chr1", ranges = IRanges(1000, 2000))
#' regions <- get_region_signal(sample.bam, regions, floor = 10)
get_region_signal <- function(sample.bam,
                                regions,
                                control.bam = NULL,
                                floor = 1,
                                read.ext = 200) {
    if (is.character(sample.bam)) {
        samp.bam <- BamFile(sample.bam)
    }

    if (is.character(control.bam)) {
        ctrl.bam <- BamFile(control.bam)
    }

    if (is.character(regions)) {
        regions <- readBed(file = regions)
    }

    samp.stat <- idxstatsBam(samp.bam)
    samp.total <- sum(samp.stat$mapped)
    samp.mmr <- samp.total / 1000000

    samp.reads <- granges(import(samp.mmr))

    if (read.ext > 0) {
        samp.reads <- extend_reads(samp.reads, downstream = read.ext)
    }

    samp.cov <- do_bedtools_coverage(a = regions, b = samp.reads, d = TRUE)

    samp.cov$coverage <- samp.cov$coverage[samp.cov$coverage > floor]
    samp.cov$coverage <- sum(samp.cov$coverage)

    samp.cov$signal <- samp.cov$coverage / samp.mmr

    regions$sample_signal <- samp.cov$signal
    metadata(regions)$sample_mmr <- samp.mmr

    ctrl.sig <- NULL
    if (!is.null(control.bam)) {
        ctrl.stat <- idxstatsBam(ctrl.bam)
        ctrl.total <- sum(ctrl.stat$mapped)
        ctrl.mmr <- ctrl.total / 1000000

        ctrl.reads <- granges(import(ctrl.bam))
        if (read.ext > 0) {
            ctrl.reads <- extend_reads(ctrl.reads, downstream = read.ext)
        }

        ctrl.cov <- do_bedtools_coverage(a = regions, b = ctrl.reads, d = TRUE)

        ctrl.cov$coverage <- ctrl.cov$coverage[ctrl.cov$coverage > floor]
        ctrl.cov$coverage <- sum(ctrl.cov$coverage)

        ctrl.cov$signal <- ctrl.cov$coverage / ctrl.mmr

        regions$control_signal <- ctrl.cov$signal
        metadata(regions)$control_mmr <- ctrl.mmr
    }

    regions
}
