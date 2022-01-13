#' Main function to get SE category
#'
#' A function to get the category and ranking of SEs in one step
#'
#' @details
#' This function is the main function used to get the category and ranking of SEs.
#' The workflow of this function is "SEfilter -> enhancerFoldchange -> SEfitspline -> SEpattern -> SEcategory"
#'
#' @param se_in pooled ROSE output file *_peaks_Gateway_SuperEnhancers.bed of all samples.
#' @param e_in pooled macs2 output file *_peaks.narrowPeak of all samples.
#' @param bl_file super-enhancer blacklist bed file download from ENCODE (ENCFF356LFX)
#' @param custom_range a vector of extra customized blacklist to ignore.
#' Format: c(chr:start-end, chr:start-end, ...).
#'
#' @param data_type quantitative file format "bam" or "bw" (default=bam)
#' @param enhancer_count_table raw count file of pooled enhancers
#' @param s1_r1_bam path of sample 1 replicate 1 bam/bw file
#' @param s1_r2_bam path of sample 1 replicate 2 bam/bw file
#' @param s2_r1_bam path of sample 2 replicate 1 bam/bw file
#' @param s2_r2_bam path of sample 2 replicate 2 bam/bw file
#' @param spline_fun spline functions ("bs", "ns", and "smooth". default=bs)
#' @param permut if you want permutation (default=TRUE)
#' @param times permutation times (default=10)
#' @param cutoff_v fold change lower and upper cutoff vector.
#' if permutation is used, it will use the value calculated by permutation. (default=c(-1.5,1.5))
#' @param s1_pair if sample 1 is paired-end (default=FALSE)
#' @param s2_pair if sample 2 is paired-end (default=FALSE)
#'
#' @return
#' A list of 7 datasets: lfc_shrink, permutation density plot, pattern plots of each SE,
#' enhancer fitted dataset, significant threshold cutoff, SE category and boxplot of SE category
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
#' # no blacklist with permutation 10 times with b-spline function
#' se_main_list <- DASE(se_in=pooled_rose,e_in=pooled_enhancer,
#' s1_r1_bam=s1_r1_path,s1_r2_bam=s1_r2_path,
#' s2_r1_bam=s2_r1_path,s2_r2_bam=s2_r2_path)
#'
#' # with blacklist and permutation 10 times with n-spline function
#' se_main_list <- DASE(se_in=pooled_rose,e_in=pooled_enhancer,bl_file= blacklist,
#' s1_r1_bam=s1_r1_path,s1_r2_bam=s1_r2_path,s2_r1_bam=s2_r1_path,s2_r2_bam=s2_r2_path,
#' spline_fun="ns")
#'
#' # no blacklist nor permutation 10 times with smooth spline function
#' se_main_list <- DASE(se_in=pooled_rose,e_in=pooled_enhancer,permut=FALSE,
#' s1_r1_bam=s1_r1_path,s1_r2_bam=s1_r2_path,s2_r1_bam=s2_r1_path,s2_r2_bam=s2_r2_path,
#' spline_fun="smooth")
#'

DASE <- function(se_in,e_in,bl_file,custom_range,
                 enhancer_count_table,data_type="bam",s1_pair=F,s2_pair=F,
                 spline_fun="bs",
                 permut=T,times=10,cutoff_v=c(-1,1),
                 s1_r1_bam,s1_r2_bam,s2_r1_bam,s2_r2_bam) {

  # Step 1: filter super-enhancer with SEfilter.R with or without SE blacklist file
  print("Step 1: merge and filter SE")
  step_1_out <- SEfilter(se_in,bl_file=bl_file,custom_range=custom_range)

  # Step 2: get enhancer fold changes within each super-enhancer with enhancerFoldchange.R
  # Merged SE Output from step 1 as SE input for step 2
  print("Step 2: calculate log2FC of constituent enhancer using Deseq2")

  merged_se_df <- step_1_out$se_merged
  if (data_type == "bam" & missing(enhancer_count_table)) {
    step_2_out <- enhancerFoldchange_bam(e_df=e_in,se_df=merged_se_df,
                                       s1_pair=s1_pair,s2_pair=s2_pair,
                                       s1_r1_bam=s1_r1_bam,s1_r2_bam=s1_r2_bam,
                                       s2_r1_bam=s2_r1_bam ,s2_r2_bam=s2_r2_bam)
  } else if (data_type=="bw" & missing(enhancer_count_table)){
    step_2_out <- enhancerFoldchange_bw(e_df=e_in,se_df=merged_se_df,
                                        s1_r1_bw=s1_r1_bam,s1_r2_bw=s1_r2_bam,
                                        s2_r1_bw=s2_r1_bam ,s2_r2_bw=s2_r2_bam)
  } else if (!missing(enhancer_count_table)) {
    step_2_out <- enhancerFoldchange_count(e_df=e_in,se_df=merged_se_df,
                                           raw_count=enhancer_count_table)
  }

  # Step 3: using weighted b-spline to fit the fold change
  if (spline_fun=="bs") {
    print("Step 3: b-spline fit log2FC")
    print(paste0("Processing total of ",length(unique(step_2_out$enhancer_deseq_result$se_merge_name))," SEs"))
    step_3_out <- SEfitspline_bs(step_2_out$enhancer_deseq_result)
  } else if(spline_fun=="ns") {
    print("Step 3: natural spline fit log2FC")
    print(paste0("Processing total of ",length(unique(step_2_out$enhancer_deseq_result$se_merge_name))," SEs"))
    step_3_out <- SEfitspline_ns(step_2_out$enhancer_deseq_result)
  } else {
    print("Step 3: smooth spline fit log2FC")
    print(paste0("Processing total of ",length(unique(step_2_out$enhancer_deseq_result$se_merge_name))," SEs"))
    step_3_out <- SEfitspline_smooth(step_2_out$enhancer_deseq_result)
  }


  # Step 4: permutation to get cutoff (normal or stringent)
  print("Step 4: permutation to get log2FC cutoff")
  e_not_in_se <- step_2_out$count_matrix_not_in_se[,c("e_merge_name",
                                                     "S1_r1","S1_r2","S2_r1","S2_r2",
                                                     "se_merge_name")]
  e_in_se <- step_2_out$enhancer_deseq_result[,c("e_merge_name",
                             "S1_r1","S1_r2","S2_r1","S2_r2",
                             "se_merge_name")]


  # chose different permut
  if (permut){
    sample_pool <- e_in_se
    step_4_out <- SEpermut(step_3_out$se_fit_df,sample_pool,times=times)
    cutoff_vector <- step_4_out$cutoff
  } else {
    step_4_out <- SEpermut(step_3_out$se_fit_df,sample_pool,permut=F)
    cutoff_vector <- cutoff_v
  }

  # step_5: segment
  print("Step 5: pattern segments process")
  step_5_out <- SEpattern(step_3_out$se_fit_df,step_3_out$spline_plot_df,
                               cutoff_vector)

  # step_6: category
  print("Step 6: final category estimate")
  step_6_out <- SEcategory(step_5_out$se_segment_percent,step_3_out$se_fit_df)


  # make final output: final category and ranking file,
  # permutation density plots
  # pattern plots of each SEs,
  # enhancer fitted fold change and counts,
  # cutoff_vector,
  # and SE segments profile to output list

  output_list <- list(lfc_shrink = step_2_out$lfc_shrink,
                     se_fit = step_3_out$se_fit_df,
                     cutoff = cutoff_vector,
                     density_plot = step_4_out$density_plot,
                     pattern_list = step_5_out$plots,
                     se_category = step_6_out$se_cat_rank,
                     boxplot =step_6_out$cat_boxplot
  )
  return(output_list)
}
