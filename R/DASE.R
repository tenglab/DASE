#' Main function to get SE category
#'
#' A function to get the category and ranking of SEs in one step
#'
#' @details
#' This function is the main function used to get the category and ranking of SEs.
#' The workflow of this function is "SEfilter -> enhancerFoldchange -> SEfitspline -> SEpattern -> SEcategory"
#'
#' @param se_in pooled SE bed file of all samples.
#' @param e_in pooled enhancer bed file of all samples.
#' @param bl_file super-enhancer blacklist bed file download from ENCODE (ENCFF356LFX)
#' @param custom_range a vector of extra customized blacklist to ignore.
#' Format: c(chr:start-end, chr:start-end, ...).
#'
#' @param data_type quantitative file format "bam" or "bw" (default=bam)
#' @param enhancer_count_table raw count file of pooled enhancers
#' @param condition_1 path of bam/bw files for condition 1 replicates, separated by ",".
#'                    Example: condition_1="path_to_bam/c1_r1.bam,path_to_bam/c1_r2.bam,...".
#' @param condition_2 path of bam/bw files for condition 1 replicates, separated by ",".
#'                    Example: condition_2="path_to_bam/c2_r1.bam,path_to_bam/c2_r2.bam,...".
#' @param c1_pair if replicates of condition 1 is paired-end, separated by ",".
#'                Order need to match with parameter condition_1.
#'                Example: c1_pair="FALSE(or F),TRUE(or T),...". (only available for bam file)
#' @param c2_pair if replicates of condition 1 is paired-end, separated by ",".
#'                Order need to match with parameter condition_2.
#'                Example: c2_pair="FALSE(or F),TRUE(or T),...". (only available for bam file)
#' @param c1_n number of replicates (samples) in condition 1 (only available for raw count table).
#' @param c2_n number of replicates (samples) in condition 2 (only available for raw count table).
#' @param spline_fun spline functions ("bs", "ns", and "smooth". default=bs)
#' @param permut if you want permutation (default=TRUE)
#' @param times permutation times (default=10)
#' @param cutoff_v fold change lower and upper cutoff vector.
#' if permutation is used, it will use the value calculated by permutation. (default=c(-1.5,1.5))
#'
#' @return
#' The output of DASE is a list with multiple data types including:
#'
#' lfc_shrink: a shrinking lfc object from _DESeq2_ for all enhancers. It can be used to creat MA plot
#' cutoff: significant threshold for fitted log2 fold changes.
#' density_plot: a density plot of permutation and original fitted log2 fold changes, if *permut=T*.
#' boxplot: a boxplot of final SE categories
#' se_category: a data frame containing final SE categories
#' pattern_list: a list containing figures for each SE pattern
#' ce_fit: a data frame containing DESeq2 output and spline-fitted log2 fold change of all constitute enhancers
#'
#' Each column in se_category represents:
#'
#' se_merge_name: name of merged SE,"chr_start_end".
#' total_width: width of merged SE (unit=k).
#' number_enhancer: number of CEs in each SE.
#' category: SE category identified by _DASE_.
#' direction: enrichment direction of SEs (none:Other or non-differential category; +: enriched in sample 2; -: enriched in sample 1; l: sample 1 shifted in 5' direction; r: sample 2 shifted in 5' direction).
#' non_mid_percent: total activity occupancy of the segments that go beyond the threshold cutoffs.
#' mean_FC: mean of log2 SE coverage fold change.
#' rank: SE category ranking based on non_mid_percent first and mean_FC for each SE category. (rank=1 means the most changed in the corresponding SE category.)
#'
#' Each column in ce_fit represents:
#'
#' e_merge_name: name of merged CEs.
#' chr: chromosome of CEs.
#' start: start posation of CEs.
#' end: end posation of CEs.
#' width: width of CEs.
#' C1_1: coverage of condition 1 replicate 1.
#' C1_2: coverage of condition 1 replicate 2.
#' C2_1: coverage of condition 2 replicate 1.
#' C2_2: coverage of condition 2 replicate 2.
#' C1_1_norm: normalized coverage of condition 1 replicate 1.
#' C1_2_norm: normalized coverage of condition 1 replicate 2.
#' C2_1_norm: normalized coverage of condition 2 replicate 1.
#' C2_2_norm: normalized coverage of condition 2 replicate 2.
#' baseMean: normalized coverage mean of all samples.
#' log2FoldChange: log2 fold change of normalized coverage.
#' lfcSE: log2 fold change standard error.
#' stat: statistic of log2 fold change.
#' pvalue: p-value of log2 fold change.
#' padj: adjusted p-value of log2 fold change.
#' baseMean_shrink: shrinkage of normalized coverage mean of all samples.
#' log2FoldChange_shrink: shrinkage of log2 fold change..
#' lfcSE_shrink: shrinkage of log2 fold change standard error.
#' pvalue_shrink: shrinkage of log2 fold change p-value
#' padj_shrink: shrinkage of log2 fold change adjusted p-value
#' se_merge_name: merged SE names
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
#' se_main_list <- DASE(se_in=pooled_se,e_in=pooled_enhancer,
#' condition_1="s1_r1_path,s1_r2_path",
#' condition_2="s2_r1_path,s2_r2_path",
#' c1_pair="F,F",c2_pair="F,F")
#'
#' # with blacklist and permutation 10 times with n-spline function
#' se_main_list <- DASE(se_in=pooled_se,e_in=pooled_enhancer,bl_file= blacklist,
#' s1_r1_bam=s1_r1_path,s1_r2_bam=s1_r2_path,s2_r1_bam=s2_r1_path,s2_r2_bam=s2_r2_path,
#' spline_fun="ns")
#'
#' # no blacklist nor permutation 10 times with smooth spline function
#' se_main_list <- DASE(se_in=pooled_se,e_in=pooled_enhancer,permut=FALSE,
#' s1_r1_bam=s1_r1_path,s1_r2_bam=s1_r2_path,s2_r1_bam=s2_r1_path,s2_r2_bam=s2_r2_path,
#' spline_fun="smooth")
#'

DASE <- function(se_in,e_in,bl_file,custom_range,
                 enhancer_count_table,data_type="bam",
                 spline_fun="bs",
                 permut=T,times=10,cutoff_v=c(-1,1),
                 condition_1,condition_2,c1_pair,c2_pair,c1_n,c2_n) {

  # Step 1: filter super-enhancer with SEfilter.R with or without SE blacklist file
  print("Step 1: merge and filter SE")
  step_1_out <- SEfilter(se_in,bl_file=bl_file,custom_range=custom_range)

  # Step 2: get enhancer fold changes within each super-enhancer with enhancerFoldchange.R
  # Merged SE Output from step 1 as SE input for step 2

  # get replicates number in each condition
  if (!missing(enhancer_count_table)) {
    c1_n <- c1_n
    c2_n <- c2_n
  } else {
    c1_bam <- unlist(strsplit(gsub(" ","",condition_1),","))
    c2_bam <- unlist(strsplit(gsub(" ","",condition_2),","))

    c1_n <- length(c1_bam)
    c2_n <- length(c2_bam)
  }

  merged_se_df <- step_1_out$se_merged
  if (data_type == "bam" & missing(enhancer_count_table)) {
    print("Step 2: calculate log2FC of constituent enhancer using Deseq2 with bam file")
    step_2_out <- enhancerFoldchange_bam(e_df=e_in,se_df=merged_se_df,
                                         condition_1=condition_1,
                                         condition_2=condition_2,
                                         c1_pair=c1_pair,
                                         c2_pair=c2_pair)

  } else if (data_type=="bw" & missing(enhancer_count_table)){
    print("Step 2: calculate log2FC of constituent enhancer using Deseq2 with bw file")
    step_2_out <- enhancerFoldchange_bw(e_df=e_in,se_df=merged_se_df,
                                        condition_1=condition_1,
                                        condition_2=condition_2)

  } else if (!missing(enhancer_count_table)) {
    print("Step 2: calculate log2FC of constituent enhancer using Deseq2 with raw enhancer count table")
    step_2_out <- enhancerFoldchange_count(e_df=e_in,se_df=merged_se_df,
                                           raw_count=enhancer_count_table,
                                           c1_n=c1_n,
                                           c2_n=c2_n)
  }

  # Step 3: using weighted b-spline to fit the fold change
  if (spline_fun=="bs") {
    print("Step 3: b-spline fit log2FC")
    print(paste0("Processing total of ",length(unique(step_2_out$enhancer_deseq_result$se_merge_name))," SEs"))
    step_3_out <- SEfitspline_bs(step_2_out$enhancer_deseq_result,c1_n=c1_n,c2_n=c2_n)
  } else if(spline_fun=="ns") {
    print("Step 3: natural spline fit log2FC")
    print(paste0("Processing total of ",length(unique(step_2_out$enhancer_deseq_result$se_merge_name))," SEs"))
    step_3_out <- SEfitspline_ns(step_2_out$enhancer_deseq_result,c1_n=c1_n,c2_n=c2_n)
  } else {
    print("Step 3: smooth spline fit log2FC")
    print(paste0("Processing total of ",length(unique(step_2_out$enhancer_deseq_result$se_merge_name))," SEs"))
    step_3_out <- SEfitspline_smooth(step_2_out$enhancer_deseq_result,c1_n=c1_n,c2_n=c2_n)
  }


  # Step 4: permutation to get cutoff
  e_in_se <- subset(step_2_out$enhancer_deseq_result,
                    select=c('e_merge_name',colnames(step_2_out$enhancer_deseq_result)[6:(5+c1_n+c2_n)],
                                                              'se_merge_name'))


  # chose different permut
  if (permut & spline_fun=="bs"){
    print("Step 4: permutation with bs-spline to get log2FC cutoff")
    sample_pool <- e_in_se
    step_4_out <- SEpermut_bs(step_3_out$se_fit_df,sample_pool,times=times,
                              c1_n=c1_n,c2_n=c2_n)
    cutoff_vector <- step_4_out$cutoff
  } else if (permut & spline_fun=="ns") {
    print("Step 4: permutation with natrual spline to get log2FC cutoff")
    sample_pool <- e_in_se
    step_4_out <- SEpermut_ns(step_3_out$se_fit_df,sample_pool,times=times,
                              c1_n=c1_n,c2_n=c2_n)
    cutoff_vector <- step_4_out$cutoff
  } else if (permut & spline_fun=="smooth") {
    print("Step 4: permutation with smooth spline to get log2FC cutoff")
    sample_pool <- e_in_se
    step_4_out <- SEpermut_smooth(step_3_out$se_fit_df,sample_pool,times=times,
                              c1_n=c1_n,c2_n=c2_n)
    cutoff_vector <- step_4_out$cutoff
  } else if (permut==FALSE){
    print("Step 4: no permutation, skipping")
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


  # change ce_fit name
  step_3_out$se_fit_df <- step_3_out$se_fit_df[,-c(26:32)]
  colnames(step_3_out$se_fit_df)[20:24] <- c("baseMean_shrink","log2FoldChange_shrink","lfcSE_shrink",
                                             "pvalue_shrink","padj_shrink")

  # make final output: final category and ranking file,
  # permutation density plots
  # pattern plots of each SEs,
  # enhancer fitted fold change and counts,
  # cutoff_vector,
  # and SE segments profile to output list

  output_list <- list(lfc_shrink = step_2_out$lfc_shrink,
                     ce_fit = step_3_out$se_fit_df,
                     cutoff = cutoff_vector,
                     density_plot = step_4_out$density_plot,
                     pattern_list = step_5_out$plots,
                     se_category = step_6_out$se_cat_rank,
                     boxplot =step_6_out$cat_boxplot
  )
  return(output_list)
}
