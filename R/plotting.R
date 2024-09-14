#' Plot enhancer ranking curve
#'
#' Plots the enhancer ranking curve based on the signal values in the regions, 
#' with visual indicators for the signal cutoff and the number of super-enhancers identified.
#'
#' @param regions A GRanges object containing `rankby_signal` and `super` columns, as well as the `threshold` value in its metadata.
#' @param style Character string specifying the plotting style.
#'   Default is "ROSE".
#' @param factor_label Optional character string for prefixing the y-axis label with a specific factor,
#'   e.g. "(H3K27ac)". 
#'   Default is `NULL`.
#'
#' @return None. This function generates a plot.
#' 
#' @author Jared Andrews, based on code from ROSE.
#' 
#' @export
#' @importFrom graphics plot abline lines text
#' @importFrom S4Vectors metadata
#'
#' @examples
#' regions <- GRanges(seqnames = "chr1", ranges = IRanges(1000, 2000))
#' regions$rankby_signal <- rnorm(length(regions))
#' regions$super <- regions$rankby_signal > 1
#' metadata(regions)$threshold <- 1
#' plot_enhancer_curve(regions)
plot_enhancer_curve <- function(regions, style = "ROSE", factor_label = NULL) {
    ylabel <- ifelse(metadata(regions)$control_subtracted, "Sample Signal - Control Signal", "Sample Signal")

    if (!is.null(factor_label)) {
        ylabel <- paste(factor_label, ylabel)
    }

    # Reversed ranking
    rev_rank <- rev(seq_len(length(regions$rankby_signal)))
    ranked_sig <- regions$rankby_signal[order(regions$rankby_signal, decreasing = TRUE)]

    plot(rev_rank, ranked_sig,
        col = "red", xlab = "Putative Enhancers",
        ylab = ylabel, pch = 19, cex = 2
    )

    abline(h = metadata(regions)$threshold, col = "grey60", lty = 2)
    abline(v = length(regions$rankby_signal) - sum(regions$super), col = "grey60", lty = 2)
    lines(rev_rank, ranked_sig, lwd = 4, col = "red")
    text(0, 0.8 * max(regions$rankby_signal), paste(
        " Cutoff used: ", metadata(regions)$threshold, "\n",
        "Super-Enhancers identified: ", sum(regions$super)
    ), pos = 4)
}
