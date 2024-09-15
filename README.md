<p align="center" width="100%">
    <img src="inst/logo/ECLIPSE_Hex.png" alt="ECLIPSE" height="330">
</p>

---

ECLIPSE (**E**nhancer **C**alling and **L**inking with **I**ntegrated **P**rofiling and **S**tructure **E**valuation) provides a performant 
implementation of the [rank ordering of super enhancers (ROSE)](http://younglab.wi.mit.edu/super_enhancer_code.html) method for identifying super enhancers.
It provides options to increase the robustness of ROSE via signal transformations prior to thresholding and additional thresholding approaches.
It also increases flexibility by exposing parameters hidden in the original implementation.
ECLIPSE additionally contains novel functionality to identify sub-structural changes within enhancer regions between groups via sliding window and binning approaches.
It also contains visualization functions to generate high-quality plots of specific loci alongside arbitrary user-provided data.

This project was conceptualized for and initially developed during the [St. Jude Children's Research Hospital KIDS24 BioHackathon](https://github.com/stjude-biohackathon) by:
- Jared Andrews (team lead)
- Alyssa Obermayer
- Nicolas Peterson
- Kelly McCastlain
- Jacqueline Myers
- Avery Bradley

It even snagged a prize for "Most Technically Impressive" project.

**This package is under active development and may break at any time.**

## Installation

Currently, the package can be installed from GitHub:

```
library(devtools)
devtools::install_github("stjude-biohackathon/ECLIPSE")
```

## Usage

Given a BAM file and a BED file of peaks, ROSE can be run with the `run_rose` function.
Optionally, a control BAM file for input or IgG from the same sample can be provided.

Alternatively, `bamFile` objects and a `GRanges` object for the peaks can be provided.

The output is a `GRanges` object containing all putative enhancers with their super enhancer designation in the `super` metadata column.

```r
library(ECLIPSE)

sample_bam <- "/path/to/bam.bam"
control_bam <- "/path/to/control.bam"
peaks <- "/path/to/peaks.bed"

putative_enhancers <- run_rose(sample.bam = sample_bam, control.bam = control_bam, peaks = peaks)

putative_enhancers
```

## Development Roadmap

- Add missing ROSE functionality.
  - TSS exclusion from stitching process.
  - Overlap of TSS of 3 or more unique genes canceling stiching for a putative enhancer. 
    - And ability to disable this process.
  - Enhancer-gene annotations (within 50 kb by default for ROSE).
    - Ability to limit to expressed genes.

## References

If you use this package in published research, please cite the original papers utilizing and describing ROSE:

[Master Transcription Factors and Mediator Establish Super-Enhancers at Key Cell Identity Genes](http://www.cell.com/abstract/S0092-8674(13)00392-9)
Warren A. Whyte, David A. Orlando, Denes Hnisz, Brian J. Abraham, Charles Y. Lin, Michael H. Kagey, Peter B. Rahl, Tong Ihn Lee and Richard A. Young
Cell 153, 307-319, April 11, 2013

[Selective Inhibition of Tumor Oncogenes by Disruption of Super-enhancers](http://www.cell.com/abstract/S0092-8674(13)00393-0)
Jakob Lovén, Heather A. Hoke, Charles Y. Lin, Ashley Lau, David A. Orlando, Christopher R. Vakoc, James E. Bradner, Tong Ihn Lee, and Richard A. Young
Cell 153, 320-334, April 11, 2013
