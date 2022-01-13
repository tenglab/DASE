#' Fit enhancer fold change and normal premut
#'
#' A function uses smooth.spline weighted by basemean to fit the enhancer fold change
#' within each SE.
#'
#' @details
#' This function will fit the enhancer fold change within each SE by using smooth.spline
#' weighted by the max basemean of two samples.
#' @param e_deseq_out_df enhancer output file from enhancerFoldchange
#'
#' @return
#' se_fit_df: enhancer fitted dataframe
#' spline_plot_df_out: spline fitted line plot dataframe
#'
#' @import DESeq2
#' @import Rsubread
#' @import apeglm
#' @import splines
#' @import ggplot2
#'
#' @export
#' @examples
#' fit_list <- SEfitspline_smooth(enhancerFoldchange_out)
#'

SEfitspline_smooth <- function(e_deseq_out_df){
  #--------------------------------------------------------------
  e_deseq_out_df$s1_mean <- (e_deseq_out_df$S1_r1_norm+e_deseq_out_df$S1_r2_norm)/2
  e_deseq_out_df$s2_mean <- (e_deseq_out_df$S2_r1_norm+e_deseq_out_df$S2_r2_norm)/2
  e_deseq_out_df$max_mean <- apply(e_deseq_out_df[,c("s1_mean","s2_mean")],1,FUN=max)

  # create original smooth.spline fit
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
      # when the number of CE less than 5, using b-spline
      if (n_top <= 4) {
        # df=2
        spline_bs_fit <- lm(log2FoldChange.1~bs(width_mid,
                                                degree=2,
                                                Boundary.knots = c(min(width_mid)-5,max(width_mid)+5)),
                            data = temp_pattern,
                            weights = max_mean)

        # save fitted.values for permutation
        temp_pattern$spline_bs <- spline_bs_fit$fitted.values

        # pattern plot
        temp_plot_df <- data.frame(se_merge_name=rep(unique(temp_pattern$se_merge_name),n),
                                   new_x=new_x$width_mid,
                                   pred_y=predict(spline_bs_fit,new_x))

      } else if (n_top >= 5) {
        # df=4
        # smooth_spline
        spline_bs_fit <- smooth.spline(x=temp_pattern$width_mid,y=temp_pattern$log2FoldChange.1,
                                       w=temp_pattern$max_mean,
                                       df=4)

        # save fitted.values for permutation
        temp_pattern$spline_bs <- spline_bs_fit$y

        # pattern plot
        temp_plot_df <- data.frame(se_merge_name=rep(unique(temp_pattern$se_merge_name),n),
                                   new_x=new_x$width_mid,
                                   pred_y=predict(spline_bs_fit,new_x)$y$width_mid)

      }




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
