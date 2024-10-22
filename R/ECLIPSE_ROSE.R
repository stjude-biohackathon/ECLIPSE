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
#' @importFrom GenomicRanges strand start end 'ranges<-'
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
#'   Coverage must be greater than this value to be counted in the signal.
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
#' @importFrom GenomicRanges granges trim
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
        sample.bam <- BamFile(sample.bam)
    }

    if (is.character(control.bam)) {
        control.bam <- BamFile(control.bam)
    }

    if (is.character(regions)) {
        regions <- readBed(file = regions)
    }

    samp_stat <- idxstatsBam(sample.bam)
    samp_total <- sum(samp_stat$mapped)
    samp_mmr <- samp_total / 1000000

    samp_reads <- granges(import(sample.bam))
    samp_reads <- trim(samp_reads)

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
        ctrl_stat <- idxstatsBam(control.bam)
        ctrl_total <- sum(ctrl_stat$mapped)
        ctrl_mmr <- ctrl_total / 1000000

        ctrl_reads <- granges(import(control.bam))
        ctrl_reads <- trim(ctrl_reads)
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
#'   Default is `TRUE`, as that is what ROSE does.
#'
#' @return A GRanges object with an added `rank_signal` column containing the
#'   computed `rank_signal` column in its metadata columns, sorted by said column.
#'   Adds a `region_rank` column as well. Adds a `control_subtracted` element to the `metadata`.
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

    metadata(regions)$control_subtracted <- FALSE
    if (!is.null(regions$control_signal)) {
        rank_sig <- rank_sig - regions$control_signal
        metadata(regions)$control_subtracted <- TRUE
    }

    if (negative.to.zero) {
        rank_sig[rank_sig < 0] <- 0
    }

    regions$rank_signal <- rank_sig
    regions <- sort(regions, decreasing = TRUE, by = ~rank_signal)
    regions$region_rank <- seq_len(NROW(regions))
    regions
}


#' Classify enhancers based on signal thresholds
#'
#' Classifies enhancers as super enhancers or typical enhancers based on a ranking signal and
#' a specified thresholding method.
#' Optionally applies user-transformations to the signal before classification, which is highly recommended
#' to ameliorate the effects of outliers on the classification.
#'
#' @param regions A GRanges object containing `rank_signal` and optionally other metadata.
#' @param transformation A function to apply to the ranking signal before threshold determination.
#'   Default is `NULL`.
#' @param drop.zeros Logical indicating whether to drop regions with zero signal.
#'   Default is `FALSE`.
#' @param thresh.method Character specifying the method to determine the signal threshold.
#'   Must be one of "ROSE", "first", "curvature", or "arbitrary".
#'   Default is "ROSE".
#' @param first.threshold Numeric value for the fraction of steepest slope when using the "first" threshold method.
#'   Higher values will result in fewer SEs called.
#'   Default is 0.5.
#' @param arbitrary.threshold Numeric value for the arbitrary threshold if the "arbitrary" method is selected.
#'   Default is 0.4, which is a reasonable setting when a cumulative proportion of signal transformation is applied.
#'
#' @return A GRanges object with a new `super` logical column indicating whether the enhancer is classified as a super enhancer.
#'   A `rankby_signal` column is added as the final signal values used for ranking (post-transformation, if applied).
#'   Any transformations applied, the thresholding method used, the threshold, and the number of dropped regions if
#'   `drop.zeros = TRUE` are added to the metadata of the GRanges object.
#'
#' @export
#'
#' @importFrom S4Vectors 'metadata<-'
#' @importFrom KneeArrower findCutoff
#'
#' @author Jared Andrews
#'
#' @examples
#' regions <- GRanges(seqnames = "chr1", ranges = IRanges(1000, 2000))
#' regions$rank_signal <- rnorm(length(regions))
#' classified_regions <- classify_enhancers(regions, thresh.method = "ROSE")
classify_enhancers <- function(regions,
                               transformation = NULL,
                               drop.zeros = FALSE,
                               thresh.method = "ROSE",
                               first.threshold = 0.5,
                               arbitrary.threshold = 0.4) {
    if (is.null(regions$rank_signal)) {
        stop("regions must contain ranking signal, run 'add_signal_rank'")
    }

    # transformation must be a function if provided
    if (!is.null(transformation) && !is.function(transformation)) {
        stop("transformation must be a function")
    }

    if (!thresh.method %in% c("ROSE", "first", "curvature", "arbitrary")) {
        stop("thresh.method must be one of 'ROSE', 'first', 'curvature', or 'arbitrary'")
    }

    metadata(regions)$threshold_method <- thresh.method

    if (drop.zeros) {
        num_no_sig_regions <- NROW(regions[regions$rank_signal == 0])
        message(paste("Dropped", num_no_sig_regions, "regions due to no signal"))
        metadata(regions)$dropped_zero_count <- num_no_sig_regions
        regions <- regions[regions$rank_signal > 0]
        regions$region_rank <- seq_len(NROW(regions))
    }

    # Keep track of whether to use transformed signal.
    use_transformed <- FALSE
    if (!is.null(transformation)) {
        message("Applying provided transformation to signal prior to determining cutoff")
        metadata(regions)$transformation <- transformation
        regions$transformed_signal <- transformation(regions$rank_signal)
        use_transformed <- TRUE
    }

    if (use_transformed) {
        regions$rankby_signal <- regions$transformed_signal
    } else {
        regions$rankby_signal <- regions$rank_signal
    }

    # Use y-axis position for threshold, i.e. the signal value rather than rank
    if (thresh.method == "ROSE") {
        cutoff_options <- calculate_cutoff(regions$rankby_signal, drawPlot = FALSE)
        cutpoint <- cutoff_options$absolute
    } else if (thresh.method == "first") {
        cutpoint <- findCutoff(
            rank(regions$rankby_signal),
            regions$rankby_signal,
            method = "first",
            frac.of.steepest.slope = first.threshold
        )
        cutpoint <- cutpoint$y
    } else if (thresh.method == "curvature") {
        cutpoint <- findCutoff(rank(regions$rankby_signal),
            regions$rankby_signal,
            method = "curvature"
        )
        cutpoint <- cutpoint$y
    } else if (thresh.method == "arbitrary") {
        cutpoint <- arbitrary.threshold
    }

    metadata(regions)$threshold <- cutpoint
    message(paste("Using", cutpoint, "as cutoff for SE classification"))

    regions$super <- regions$rankby_signal > cutpoint
    message(paste(sum(regions$super), "super enhancers called"))

    regions
}


#' Run ROSE (Rank Ordering of Super-Enhancers)
#'
#' This function performs the ROSE for identifying super-enhancers by stitching
#' together peaks, calculating the signal in the regions, ranking the regions by signal,
#' and classifying them as super-enhancers.
#'
#' @details
#' This function allows for near identical functionality to the original ROSE software,
#' but also provides a number of additional options to minimize the impact of signal outliers
#' on the classification threshold.
#' In particular, the `transformation` argument allows for
#' the application of a function to the ranking signal before threshold determination.
#' This is highly recommended to prevent outliers from skewing the threshold.
#'
#' The `thresh.method` argument allows for the selection of the method to determine the threshold.
#' - The "ROSE" method is the default and uses a sliding diagonal line to determine the cutoff.
#' - The "first" method uses a first derivative to determine the point where the slope is a given fraction of the maximum.
#'   The `first.threshold` argument controls this fraction.
#' - The "curvature" method finds the point at which the circle tanget to the curve has the smallest radius.
#' - The "arbitrary" method allows for the user to specify a fixed threshold, useful for transformations
#'   that result in a consistent curve shape with a known maximum, like cumulative proportion of signal.
#'
#' @param sample.bam A character string or `BamFile` object representing the sample BAM file.
#' @param peaks A character string or `GRanges` object representing the peaks.
#' @param control.bam A character string or `BamFile` object representing the control BAM file.
#'   Default is `NULL`.
#' @param stitch.distance Numeric value for the distance within which peaks are stitched together.
#'   Default is 12500.
#' @param exclude.tss.distance Numeric value for distance from TSS to exclude peaks prior to stitching.
#'   Default is 0, meaning peaks overlapping TSS will not be excluded from stitching.
#' @param txdb Not used. Functionality not implemented.
#'   Default is `NULL`.
#' @param org.db Not used. Functionality not implemented.
#'   Default is `NULL`.
#' @param drop.y Logical indicating whether to drop peaks on chromosome Y.
#'   Default is `TRUE`.
#' @param max.genes.overlap Not used. Functionality not implemented.
#'   Default is 2.
#' @param negative.to.zero Logical indicating whether to set negative ranking signals to zero.
#'   Default is `TRUE`.
#' @param thresh.method Character string specifying the method to determine the signal threshold.
#'   Must be one of "ROSE", "first", "curvature", or "arbitrary".
#' Default is "ROSE".
#' @param transformation A function to apply to the ranking signal before threshold determination.
#'   Default is `NULL`.
#' @param floor Numeric value representing the minimum coverage threshold to count.
#'  Default is 1.
#' @param read.ext Numeric value for extending reads downstream.
#'   Default is 200.
#' @param drop.zeros Logical indicating whether to drop regions with zero signal.
#'   Default is `FALSE`.
#' @param first.threshold Numeric value for the fraction of steepest slope when using the "first" threshold method.
#'   Default is 0.5.
#' @param arbitrary.threshold Numeric value for the arbitrary threshold if the "arbitrary" method is selected.
#'   Default is 0.4.
#'
#' @return A `GRanges` object containing the classified super-enhancers and associated metadata.
#'
#' @author Jared Andrews
#'
#' @importFrom Rsamtools BamFile indexBam
#' @importFrom genomation readBed
#' @importFrom GenomicRanges reduce seqnames trim
#'
#' @export
#'
#' @examples
#' sample.bam <- "path/to/sample.bam"
#' peaks <- "path/to/peaks.bed"
#' result <- run_rose(sample.bam, peaks)
run_rose <- function(
    sample.bam,
    peaks,
    control.bam = NULL,
    stitch.distance = 12500,
    exclude.tss.distance = 0,
    txdb = NULL, # Not used/functionality not implemented
    org.db = NULL, # Not used/functionality not implemented
    drop.y = TRUE,
    max.genes.overlap = 2, # Not used/functionality not implemented
    negative.to.zero = TRUE,
    thresh.method = "ROSE",
    transformation = NULL,
    floor = 1,
    read.ext = 200,
    drop.zeros = FALSE,
    first.threshold = 0.5,
    arbitrary.threshold = 0.4) {

    if (is.character(sample.bam)) {
        sample.bam <- BamFile(sample.bam)
    }
    message(paste0("Sample BAM file: ", sub(".*/(.*\\.bam)$", "\\1", sample.bam$path)))

    if(length(sample.bam$index) == 0 || !file.exists(sample.bam$index)) {
        message("Sample BAM index not found. Generating an index.")
        sample.bam$index <- unname(indexBam(sample.bam))
    }

    if (!is.null(control.bam)) {
        if(is.character(control.bam)) {
        control.bam <- BamFile(control.bam)
        }
        message(paste0("Control BAM file: ", sub(".*/(.*\\.bam)$", "\\1", control.bam$path)))

        if(length(control.bam$index) == 0 || !file.exists(control.bam$index)) {
            message("Control BAM index not found. Generating an index.")
            control.bam$index <- unname(indexBam(control.bam))
        }
    }

    message("Reading peaks")
    if (is.character(peaks)) {
        peaks <- readBed(peaks)
        peaks <- trim(peaks)
    }

    peaks_stitched <- reduce(peaks, min.gapwidth = stitch.distance)

    # Drop chrY as ROSE does
    if (drop.y) {
        peaks_stitched <- peaks_stitched[seqnames(peaks_stitched) != "chrY"]
    }

    message("Calculating normalized region signal")
    regions <- add_region_signal(sample.bam, peaks_stitched, control.bam = control.bam, floor = floor, read.ext = read.ext)

    message("Ranking regions")
    regions <- add_signal_rank(regions, negative.to.zero = negative.to.zero)

    message("Classifying enhancers")
    regions <- classify_enhancers(regions,
        transformation = transformation, drop.zeros = drop.zeros,
        thresh.method = thresh.method, first.threshold = first.threshold, arbitrary.threshold = arbitrary.threshold
    )

    regions
}
