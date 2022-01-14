#' Fit enhancer fold change and normal premut
#'
#' A function uses b-spline weighted by basemean to fit the enhancer fold change
#' within each SE.
#'
#' @details
#' This function will fit the enhancer fold change within each SE by using b-spline
#' weighted by the max basemean of two samples.
#' @param e_deseq_out_df enhancer output file from enhancerFoldchange
#'
#' @return
#' se_fit_df: enhancer fitted dataframe
#' spline_plot_df_out: spline fitted line plot dataframe
#'
#' Each column in se_fit_df represents:
#'
#' e_merge_name: name of merged CEs.
#' chr: chromosome of CEs.
#' start: start posation of CEs.
#' end: end posation of CEs.
#' width: width of CEs.
#' S1_r1: coverage of sample 1 replicate 1.
#' S1_r2: coverage of sample 1 replicate 2.
#' S2_r1: coverage of sample 2 replicate 1.
#' S2_r2: coverage of sample 2 replicate 2.
#' S1_r1_norm: normalized coverage of sample 1 replicate 1.
#' S1_r2_norm: normalized coverage of sample 1 replicate 2.
#' S2_r1_norm: normalized coverage of sample 2 replicate 1.
#' S2_r2_norm: normalized coverage of sample 2 replicate 2.
#' baseMean: normalized coverage mean of all samples.
#' log2FoldChange: log2 fold change of normalized coverage.
#' lfcSE: log2 fold change standard error.
#' stat: statistic of log2 fold change.
#' pvalue: p-value of log2 fold change.
#' padj: adjusted p-value of log2 fold change.
#' baseMean_shrink: shrinkage of normalized coverage mean of all samples based on the methods of _DESeq2_.
#' log2FoldChange_shrink: shrinkage of log2 fold change..
#' lfcSE_shrink: shrinkage of log2 fold change standard error.
#' pvalue_shrink: shrinkage of log2 fold change p-value
#' padj_shrink: shrinkage of log2 fold change adjusted p-value
#' se_merge_name: merged SE names
#'
#' @import DESeq2
#' @import Rsubread
#' @import apeglm
#' @import splines
#' @import ggplot2
#'
#' @export
#' @examples
#' fit_list <- SEfitspline_bs(enhancerFoldchange_out)
#'

SEfitspline_bs <- function(e_deseq_out_df){
  #--------------------------------------------------------------
  e_deseq_out_df$s1_mean <- (e_deseq_out_df$S1_r1_norm+e_deseq_out_df$S1_r2_norm)/2
  e_deseq_out_df$s2_mean <- (e_deseq_out_df$S2_r1_norm+e_deseq_out_df$S2_r2_norm)/2
  e_deseq_out_df$max_mean <- apply(e_deseq_out_df[,c("s1_mean","s2_mean")],1,FUN=max)

  # create original b-spline fit
  fit_out <- data.frame()
  spline_plot_df_out <- data.frame()

  merged_super_name <- unique(e_deseq_out_df$se_merge_name)
  for (i in c(1:length(merged_super_name))) {

    # print step information
    if(i %% 500==0) {
      # Print on the screen some message
      print(paste0("SE: ",i))
    }

    temp_pattern <- e_deseq_out_df[which(e_deseq_out_df$se_merge_name == merged_super_name[i]),]
    temp_pattern$width_mid <- (temp_pattern$start+temp_pattern$width/2)/1000
    temp_pattern$percent <- temp_pattern$max_mean/sum(temp_pattern$max_mean)*100
    temp_pattern <- temp_pattern[order(-temp_pattern$percent)]
    temp_pattern$cumsum <- cumsum(temp_pattern$percent)

    # check if only 1 enhancer
    if (nrow(temp_pattern) == 1) {
      spline_bs_fit <- lm(log2FoldChange.1~bs(width_mid,
                                              degree=2),
                          data = temp_pattern,
                          weights = max_mean)

      temp_pattern$spline_bs <- spline_bs_fit$fitted.values

      # pattern plot
      temp_plot_df <- data.frame(se_merge_name=temp_pattern$se_merge_name,
                                 new_x=temp_pattern$width_mid,
                                 pred_y=temp_pattern$log2FoldChange.1)

    } else {
      # get df based on number of enhancer whose cumsum is > 95%
      n_top <- min(which(temp_pattern$cumsum > 95))

      # create new_x
      n<-20
      new_x <- data.frame(width_mid=seq(min(temp_pattern$width_mid),
                                        max(temp_pattern$width_mid),
                                        length.out=n))

      # different df based on n_top
      if (n_top <= 4) {
        # df=2
        spline_bs_fit <- lm(log2FoldChange.1~bs(width_mid,
                                                degree=2,
                                                Boundary.knots = c(min(width_mid)-5,max(width_mid)+5)),
                            data = temp_pattern,
                            weights = max_mean)

      } else if (n_top >= 6) {
        # df=4
        spline_bs_fit <- lm(log2FoldChange.1~bs(width_mid,
                                                degree=4,
                                                Boundary.knots = c(min(width_mid)-5,max(width_mid)+5)),
                            data = temp_pattern,
                            weights = max_mean)

      } else {
        # df=3
        spline_bs_fit <- lm(log2FoldChange.1~bs(width_mid,
                                                degree=3,
                                                Boundary.knots = c(min(width_mid)-5,max(width_mid)+5)),
                            data = temp_pattern,
                            weights = max_mean)

      }

      # save fitted.values for permutation
      temp_pattern$spline_bs <- spline_bs_fit$fitted.values

      # pattern plot
      temp_plot_df <- data.frame(se_merge_name=rep(unique(temp_pattern$se_merge_name),n),
                                   new_x=new_x$width_mid,
                                   pred_y=predict(spline_bs_fit,new_x))

    }
    # create output dataframe with fitted value
    fit_out <- rbind(fit_out,temp_pattern)
    spline_plot_df_out <- rbind(spline_plot_df_out,temp_plot_df)
  }

  # creat final output list
  out <- list()
  out$se_fit_df <- fit_out
  out$spline_plot_df <- spline_plot_df_out
  return(out)
}
