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
#' @param s1_r1_in_bam path of sample 1 replicate 1 bam file
#' @param s1_r2_in_bam path of sample 1 replicate 2 bam file
#' @param s2_r1_in_bam path of sample 2 replicate 1 bam file
#' @param s2_r2_in_bam path of sample 2 replicate 2 bam file
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

SEprofile_w_input <- function(se_in,e_df,bl_file=FALSE,has_bl_file=FALSE,
                      s1_pair=FALSE,s2_pair=FALSE, bw=F, custom_range=F,
                      permut=TRUE, permut_type="normal",times=10,cutoff_v=c(-1,1),
                      s1_r1_bam,s1_r2_bam,s2_r1_bam,s2_r2_bam,
                      s1_r1_in_bam,s1_r2_in_bam,s2_r1_in_bam,s2_r2_in_bam) {

  # Step 1: filter super-enhancer with SEfilter.R with or without SE blacklist file
  print("Step 1: merge and filter SE")
  step_1_out <- SEfilter(se_in,bl_file,has_bl_file,custom_range)

  # Step 2: get enhancer fold changes within each super-enhancer with enhancerFoldchange.R
  # Merged SE Output from step 1 as SE input for step 2
  print("Step 2: calculate log2FC of constituent enhancer using Deseq2")

  merged_se_df <- step_1_out$se_merged
  print(paste0("Processing total of ",nrow(merged_se_df)," SEs"))
  if (bw == F) {
    step_2_out <- enhancerFoldchange_w_input(e_df=e_df,se_df=merged_se_df,
                                       s1_pair=s1_pair,s2_pair=s2_pair,
                                       s1_r1_bam=s1_r1_bam,s1_r2_bam=s1_r2_bam,
                                       s2_r1_bam=s2_r1_bam ,s2_r2_bam=s2_r2_bam,
                                       s1_r1_in_bam = s1_r1_in_bam_path,
                                       s1_r2_in_bam = s1_r2_in_bam_path,
                                       s2_r1_in_bam = s2_r1_in_bam_path,
                                       s2_r2_in_bam = s2_r2_in_bam_path)
  } else {
    step_2_out <- enhancerFoldchange_bw(e_df=e_df,se_df=merged_se_df,
                                        s1_pair=s1_pair,s2_pair=s2_pair,
                                        s1_r1_bam=s1_r1_bam,s1_r2_bam=s1_r2_bam,
                                        s2_r1_bam=s2_r1_bam ,s2_r2_bam=s2_r2_bam)
  }

  # Step 3: using weighted bs-spline to fit the fold change
  print("Step 3: bs-spline fit log2FC")
  step_3_out <- SEfitspline_new(step_2_out$enhancer_deseq_result)
  # Step 4: permutation to get cutoff (normal or stringent)
  print("Step 4: permutation to get log2FC cutoff")
  e_not_in_se <- step_2_out$count_matrix_not_in_se[,c("e_merge_name",
                                                     "S1_r1","S1_r2","S2_r1","S2_r2",
                                                     "se_merge_name")]
  e_in_se <- step_2_out$enhancer_deseq_result[,c("e_merge_name",
                             "S1_r1","S1_r2","S2_r1","S2_r2",
                             "se_merge_name")]


  # chose different permut_type
  if (permut_type == "normal") {
    if (permut){
      sample_pool <- e_in_se
      step_4_out <- permut_spline_new(step_3_out$se_fit_df,sample_pool,times=times)
      cutoff_vector <- step_4_out$cutoff
    } else {
      step_4_out <- permut_spline_new(step_3_out$se_fit_df,sample_pool,permut=F)
      cutoff_vector <- cutoff_v
    }

  } else if (permut_type == "stringent") {

    if (permut){
      sample_pool <- rbind(e_not_in_se,e_in_se)
      step_4_out <- permut_spline_new(step_3_out$se_fit_df,sample_pool,times=times)
      cutoff_vector <- step_4_out$cutoff
    } else {
      step_4_out <- permut_spline_new(step_3_out$se_fit_df,sample_pool,permut=F)
      cutoff_vector <- cutoff_v
    }
  }


  # step_5: segment
  print("Step 5: pattern segments process")
  step_5_out <- SEpattern_demo(step_3_out$se_fit_df,step_3_out$spline_plot_df,
                               cutoff_vector)

  # step_6: category
  print("Step 6: final category estimate")
  step_6_out <- SEcategory_demo_2(step_5_out$se_segment_percent,step_3_out$se_fit_df)


  # make final output: final category and ranking file,
  # permutation density plots
  # pattern plots of each SEs,
  # enhancer fitted fold change and counts,
  # cutoff_vector,
  # and SE segments profile to output list

  output_list <- list(se_meta = step_1_out$se_merged_meta,
                     se_deseq_out = step_2_out$enhancer_deseq_result,
                     e_not_in_se = step_2_out$count_matrix_not_in_se[,c("e_merge_name",
                                                                        "S1_r1","S1_r2","S2_r1","S2_r2",
                                                                        "se_merge_name")],
                     e_in_se = step_2_out$enhancer_deseq_result[,c("e_merge_name",
                                                                   "S1_r1","S1_r2","S2_r1","S2_r2",
                                                                   "se_merge_name")],
                     lfc_shrink = step_2_out$lfc_shrink,
                     se_fit = step_3_out$se_fit_df,
                     pattern_plot_df = step_3_out$spline_plot_df,
                     cutoff = cutoff_vector,
                     density_plot = step_4_out$density_plot,
                     pattern_list = step_5_out$plots,
                     se_seg = step_5_out$se_segment_percent,
                     se_category = step_6_out$se_cat_rank,
                     boxplot =step_6_out$cat_boxplot
  )
  return(output_list)
}
