% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enhancerFoldchange_bw.R
\name{enhancerFoldchange_bw}
\alias{enhancerFoldchange_bw}
\title{Calculate the fold change of enhancers}
\usage{
enhancerFoldchange_bw(e_df, se_df, condition_1, condition_2)
}
\arguments{
\item{e_df}{pooled macs2 output file *_peaks.narrowPeak of all samples.}

\item{se_df}{merged SE metadata from SEfilter}

\item{condition_1}{path of bam files for condition 1 replicates, separated by ",".
Example: condition_1="path_to_bw/c1_r1.bw,path_to_bw/c1_r2.bw,...".}

\item{condition_2}{path of bam files for condition 1 replicates, separated by ",".
Example: condition_2="path_to_bw/c2_r1.bw,path_to_bw/c2_r2.bw,...".}
}
\value{
A list of 2 datasets: enhancer foldchang results dataset and sizefactor used in DESeq2
}
\description{
A function to calculate the fold change of enhancers base on two conditions.
}
\details{
This function will extract the enhancers which in the range of SEs.
Then enhancers are merged if there gaps are within 500bps.
Feature count is used if provide bam file
rtracklayer is used if privide bw file
DESeq2 are used to calculate the fold change of those enhancers
}
\examples{
foldchange_list <- enhancerFoldchange(pool_enhancer_df,se_meta,
                       condition_1="path_to_bw/c1_r1.bw,path_to_bam/c1_r2.bw,...",
                       condition_2="path_to_bw/c2_r1.bw,path_to_bam/c2_r2.bw,...")

}
