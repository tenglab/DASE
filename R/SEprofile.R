#' Get SE category by one step
#'
#' A function to get the category and ranking of SEs in one step
#'
#' @details
#' This function is the main function used to get the category and ranking of SEs.
#' The workflow of this function is "SEfilter -> enhancerFoldchange -> SEfitspline -> SEpattern -> SEcategory"
#'
#' @param se_in pooled ROSE output file *_peaks_Gateway_SuperEnhancers.bed of all samples.
#' Header needs to be "CHROM,START,STOP,REGION_ID,Signal,STRAND"
#' @param e_df pooled macs2 output file *_peaks.narrowPeak of all samples.
#' Header needs to be "chr, start,end,name,score,strand,signalValue,pValue,qValue,peak"
#' @param bl_file super-enhancer blacklist bed file download from ENCODE (ENCFF356LFX) (default=FALSE)
#' @param has_bl_file if there is a blacklist file (default=FALSE)
#' @param custom_range a vector of extra customized blacklist to ignore.
#' Format: c(chr:start-end, chr:start-end, ...) (default=FALSE).
#'
#' @param bw quantitative file are bigwig files (default=FALSE)
#' @param s1_r1_bam path of sample 1 replicate 1 bam/bw file
#' @param s1_r2_bam path of sample 1 replicate 2 bam/bw file
#' @param s2_r1_bam path of sample 2 replicate 1 bam/bw file
#' @param s2_r2_bam path of sample 2 replicate 2 bam/bw file
#' @param permut if you want permutation (default=TRUE)
#' @param permut_type normal or stringent (default=normal)
#' @param times permutation times (default=10)
#' @param cutoff_v fold change lower and upper cutoff vector.
#' if permutation is used, it will use the value calculated by permutation. (default=c(-1.5,1.5))
#' @param s1_pair if sample 1 is paired-end (default=FALSE)
#' @param s2_pair if sample 2 is paired-end (default=FALSE)
#'
#' @return
#' A list of 6 datasets: permutation density plot, pattern plots of each SE,
#' enhancer fitted dataset, fold change cutoff vector,
#' SE segment percentage dataset, and SE category./ranking dataset
#'
#' @import rtracklayer
#' @import data.table
#' @import ggplot2
#' @import DESeq2
#' @import Rsubread
#' @import apeglm
#' @import splines
#'
#' @export
#' @examples
#' # no blacklist with permutation 10 times
#' se_main_list <- SEprofile(se_in=pooled_rose,e_df=pooled_enhancer,
#' s1_r1_bam=s1_r1_path,s1_r2_bam=s1_r2_path,
#' s2_r1_bam=s2_r1_path,s2_r2_bam=s2_r2_path)
#'
#' # with blacklist and permutation 10 times
#' se_main_list <- SEprofile(se_in=pooled_rose,e_df=pooled_enhancer,bl_file= blacklist,
#' has_bl_file=TRUE, s1_r1_bam=s1_r1_path,s1_r2_bam=s1_r2_path,s2_r1_bam=s2_r1_path,s2_r2_bam=s2_r2_path)
#'
#' # no blacklist nor permutation 10 times
#' se_main_list <- SEprofile(se_in=pooled_rose,e_df=pooled_enhancer,permut=FALSE,
#' s1_r1_bam=s1_r1_path,s1_r2_bam=s1_r2_path,s2_r1_bam=s2_r1_path,s2_r2_bam=s2_r2_path)
#'

SEprofile <- function(se_in,e_df,bl_file=FALSE,has_bl_file=FALSE,
                      s1_pair=FALSE,s2_pair=FALSE, bw=F, custom_range=F,
                      permut=TRUE, permut_type="normal",times=10,cutoff_v=c(-1.5,1.5),
                      s1_r1_bam,s1_r2_bam,s2_r1_bam,s2_r2_bam) {

  # Step 1: filter super-enhancer with SEfilter.R with or without SE blacklist file
  step_1_out <- SEfilter(se_in,bl_file,has_bl_file,custom_range)

  # Step 2: get enhancer fold changes within each super-enhancer with enhancerFoldchange.R
  # Merged SE Output from step 1 as SE input for step 2
  merged_se_df <- step_1_out$se_merged

  if (bw == F) {
    if (permut_type == "normal") {
      step_2_out <- enhancerFoldchange(e_df=e_df,se_df=merged_se_df,
                                       s1_pair=s1_pair,s2_pair=s2_pair,
                                       s1_r1_bam=s1_r1_bam,s1_r2_bam=s1_r2_bam,
                                       s2_r1_bam=s2_r1_bam ,s2_r2_bam=s2_r2_bam)
    } else {
      step_2_out <- enhancerFoldchange(e_df=e_df,se_df=merged_se_df,
                                       s1_pair=s1_pair,s2_pair=s2_pair,
                                       permut_type = "stringent",
                                       s1_r1_bam=s1_r1_bam,s1_r2_bam=s1_r2_bam,
                                       s2_r1_bam=s2_r1_bam ,s2_r2_bam=s2_r2_bam)
    }
  } else {
    if (permut_type == "normal") {
      step_2_out <- enhancerFoldchange_bw(e_df,merged_se_df,
                                          s1_r1_bam,s1_r2_bam,s2_r1_bam,s2_r2_bam)
    } else {
      step_2_out <- enhancerFoldchange_bw(e_df,merged_se_df,
                                          permut_type = "stringent",
                                          s1_r1_bam,s1_r2_bam,s2_r1_bam,s2_r2_bam)
    }
  }

  # Step 3: using weighted bs-spline to fit the fold change
  # and use permutation to get fold change significant cutoff
  e_deseq_out <- step_2_out$enhancer_deseq_result
  not_in_se_count_matrix <- test_2_out$count_matrix_not_in_se

  # chose different permut_type
  if (permut_type == "normal") {
    if (permut == TRUE){
      step_3_out <- SEfitspline(e_deseq_out,merged_se_df)
      cutoff_vector <- step_3_out$cutoff
    } else {
      step_3_out <- SEfitspline(e_deseq_out,merged_se_df,permut = F)
      cutoff_vector <- cutoff_v
    }

    se_fit_df <- step_3_out$se_fit_df
  } else if (permut_type == "stringent") {

    if (permut == TRUE){
      step_3_out <- stringent_permut(not_in_se_count_matrix,e_deseq_out_df,merged_se_df)
      cutoff_vector <- step_3_out$cutoff
    } else {
      step_3_out <- stringent_permut(not_in_se_count_matrix,e_deseq_out_df,merged_se_df,permut = F)
      cutoff_vector <- cutoff_v
    }

    se_fit_df <- step_3_out$se_fit_df

  }


  # Step 4: using cutoff to get the segment pattern of each SE
  step_4_out <- SEpattern(se_fit_df,cutoff_vector)
  se_seg_df <- step_4_out$se_segment_percent

  # Step 5: get categories of each SE and ranking by category groups
  step_5_out <- SEcategory(se_seg_df,se_fit_df)

  # make final output: final category and ranking file,
  # permutation density plots
  # pattern plots of each SEs,
  # enhancer fitted fold change and counts,
  # cutoff_vector,
  # and SE segments profile to output list

  out <- list()
  out$lfc_shrink <- step_2_out$lfc_shrink
  out$se_meta <- step_1_out$se_merged
  out$permut_plot <- step_3_out$density_plot
  out$e_fit_fc <- step_3_out$se_fit_df
  out$fc_cutoff <- cutoff_vector
  out$pattern_plot <- step_4_out$plots
  out$se_segments <- step_4_out$se_segment_percent
  out$cate_rank <- step_5_out$se_cat_rank
  out$box_plot <- step_5_out$cat_boxplot
  return(out)
}
