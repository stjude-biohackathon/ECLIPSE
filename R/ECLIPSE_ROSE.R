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


#' Get region signal from BAM files and add to GRanges
#'
#' Calculates the libary size-normalized coverage of reads within specified regions from a
#' signal BAM file with optional subtraction of a control BAM file.
#' Allows read extension and additional processing like setting negative coverage to zero.
#'
#' @details The defaults of this function are set to so as to match the functionality
#' of ROSE as closely as possible.
#'
#' To detail the process:
#' - The total number of reads in the signal BAM file is calculated.
#' - The signal reads are extended downstream by a specified number of bases (200 bp by default).
#' - The coverage for each basepair in the regions of interest are calculated.
#'     ROSE does this manually by calling samtools for each region, which is slow.
#' - Basepairs with coverage below a specified threshold (`floor`, 1 by default) are removed.
#' - For each region, the coverage is summed and divided by the total number of reads to get the signal.
#'
#' @param sample.bam Character or BamFile object representing the signal BAM file.
#' @param regions Character or GRanges object representing genomic regions of interest.
#' @param control.bam Character or BamFile object representing the control BAM file.
#'   Default is `NULL`.
#' @param floor Numeric value for the minimum coverage threshold to count for region.
#'   Default is 1.
#' @param read.ext Numeric value for extending reads downstream. Default is 200.
#'
#' @return A GRanges object for `regions` with additional columns for sample and control (if provided) signal.
#'   The metadata of the object will also contain the scaling factor for library size normalization for the `sample.bam`
#'   and `control.bam` (if provided) in the `sample_mmr` and `control_mmr` elements.
#'
#' @export
#' 
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
#' regions <- add_region_signal(sample.bam, regions, floor = 10)
add_region_signal <- function(sample.bam,
                              regions,
                              control.bam = NULL,
                              floor = 1,
                              read.ext = 200) {
    if (is.character(sample.bam)) {
        samp_bam <- BamFile(sample.bam)
    }

    if (is.character(control.bam)) {
        ctrl_bam <- BamFile(control.bam)
    }

    if (is.character(regions)) {
        regions <- readBed(file = regions)
    }

    samp_stat <- idxstatsBam(samp_bam)
    samp_total <- sum(samp_stat$mapped)
    samp_mmr <- samp_total / 1000000

    samp_reads <- granges(import(samp_mmr))

    if (read.ext > 0) {
        samp_reads <- extend_reads(samp_reads, downstream = read.ext)
    }

    samp_cov <- do_bedtools_coverage(a = regions, b = samp_reads, d = TRUE)

    samp_cov$coverage <- samp_cov$coverage[samp_cov$coverage > floor]
    samp_cov$coverage <- sum(samp_cov$coverage)

    samp_cov$signal <- samp_cov$coverage / samp_mmr

    regions$sample_signal <- samp_cov$signal
    metadata(regions)$sample_mmr <- samp_mmr

    ctrl_sig <- NULL
    if (!is.null(control.bam)) {
        ctrl_stat <- idxstatsBam(ctrl_bam)
        ctrl_total <- sum(ctrl_stat$mapped)
        ctrl_mmr <- ctrl_total / 1000000

        ctrl_reads <- granges(import(ctrl_bam))
        if (read.ext > 0) {
            ctrl_reads <- extend_reads(ctrl_reads, downstream = read.ext)
        }

        ctrl_cov <- do_bedtools_coverage(a = regions, b = ctrl_reads, d = TRUE)

        ctrl_cov$coverage <- ctrl_cov$coverage[ctrl_cov$coverage > floor]
        ctrl_cov$coverage <- sum(ctrl_cov$coverage)

        ctrl_cov$signal <- ctrl_cov$coverage / ctrl_mmr

        regions$control_signal <- ctrl_cov$signal
        metadata(regions)$control_mmr <- ctrl_mmr
    }

    regions
}


#' Get ranking signal for regions and add to GRanges object
#'
#' Computes the ranking signal for the given genomic regions by
#' optionally subtracting control signal and setting negative values to zero.
#'
#' @param regions A GRanges object containing the `sample_signal` and optionally `control_signal` in metadata columns.
#' @param negative.to.zero Logical indicating whether to set negative values in the ranking signal to zero. 
#'   Default is `TRUE`.
#'
#' @return A GRanges object with an added `rank_signal` column containing the
#'   computed `rank_signal` column in its metadata columns, sorted by said column.
#'   Adds a `region_rank` column as well.
#'
#' @export
#' 
#' @importFrom GenomicRanges sort
#'
#' @author Jared Andrews, Jacqueline Myers
#'
#' @examples
#' regions <- GRanges(seqnames = "chr1", ranges = IRanges(1000, 2000))
#' regions$sample_signal <- rnorm(length(regions))
#' regions$control_signal <- rnorm(length(regions))
#' ranked_regions <- get_ranking_signal(regions)
add_signal_rank <- function(regions, negative.to.zero = TRUE) {
    if (is.null(regions$sample_signal)) {
        stop("regions must contain signal, run 'get_region_signal'")
    }

    rank_sig <- regions$sample_signal

    if (!is.null(regions$control_signal)) {
        rank_sig <- rank_sig - regions$control_signal
    }

    if (negative.to.zero) {
        rank_sig[rank_sig < 0] <- 0
    }

    regions$rank_signal <- rank_sig
    regions <- sort(regions, decreasing = TRUE, by = ~ rank_signal)
    regions$region_rank <- seq_len(NROW(regions))
    regions
}
