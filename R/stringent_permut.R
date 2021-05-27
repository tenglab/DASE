#' Fit enhancer fold change and stringent premut with all enhancers
#'
#' A function uses bs-spline weighted by basemean to fit the enhancer fold change
#' within each SE.
#'
#' @details
#' This function will fit the enhancer fold change within each SE by using bs-spline
#' weighted by basemean. A stringent permutation of counts of all enhancers is used to decide the
#' cutoff of significant log2FC (for samples are too similar). Cutoff was calculate by slope equal to range(y)/range(x).
#' The max value in range(-2,-1.5) and min value in range(1.5,2) are set to be lower and upper cutoff.
#'
#' @param not_in_se_count_matrix enhancer count matrix not in SE range
#' (header "e_merge_name","S1_r1","S1_r2","S2_r1","S2_r2","se_merge_name","chr","start","end")
#' @param e_deseq_out_df enhancer output file from enhancerFoldchange
#' @param se_df merged SE metadata from SEfilter
#' @param times permutation times (default=10)
#' @param permut if you want permutation (default=TRUE)
#'
#' @return
#' A list of 4 datasets: enhancer fitted dataset, permutation dataset, cutoff vector,
#' and a density plot of permutation and original fitted fold change.
#'
#' @import data.table
#' @import DESeq2
#' @import Rsubread
#' @import apeglm
#' @import ggplot2
#'
#' @export
#' @examples
#' #' # permutation 10 times
#' fit_list <- stringent_permut(pool_enhancer_df,enhancerFoldchange_out,se_meta,
#' s1_r1_bam,s1_r2_bam,s2_r1_bam,s2_r2_bam)
#'
#' # permutation 5 times
#' fit_list <- stringent_permut(pool_enhancer_df,enhancerFoldchange_out,se_meta,times=5,
#' s1_r1_bam,s1_r2_bam,s2_r1_bam,s2_r2_bam)
#'
#' # no permutation
#' fit_list <- stringent_permut(pool_enhancer_df,enhancerFoldchange_out,se_meta,permut=FALSE,
#' s1_r1_bam,s1_r2_bam,s2_r1_bam,s2_r2_bam)
#'

stringent_permut <- function(not_in_se_count_matrix,e_deseq_out_df,se_df,
                             permut=T,times=10) {

  #--------------------------------------------------------------
  # create original bs-spline fit
  fit_out <- data.frame()
  merged_super_name <- se_df$se_merge_name
  for (i in c(1:length(merged_super_name))) {
    temp_pattern <- e_deseq_out_df[which(e_deseq_out_df$se_merge_name == merged_super_name[i]),]
    temp_pattern$width_mid <- (temp_pattern$start+temp_pattern$width/2)/1000
    # splines bs fit
    spline_bs_fit <- lm(formula = log2FoldChange.1 ~ bs(width_mid),data = temp_pattern,
                        weights = baseMean.1)

    temp_pattern$spline_bs <- spline_bs_fit$fitted.values

    # create output dataframe with fitted value
    fit_out <- rbind(fit_out,temp_pattern)
  }

  #--------------------------------------------------------------
  # make in se count matrix
  in_se_count_temp <- e_deseq_out_df[,c(1,6:9,25,2:4)]
  in_se_count_matrix <- data.frame(in_se_count_temp, row.names=in_se_count_temp$e_merge_name)

  # create total enhancer df
  all_enhancer_count_matrix <- rbind(in_se_count_matrix,not_in_se_count_matrix)

  #-----------------------------------------------------------------------------
  # shuffle counts of all enhancers n times (default 10), got new fold change within SEs,
  # and get new bs-spline fit
  merged_super_name <- se_df$se_merge_name

  in_permut_out <- data.frame()
  all_permut_out <- data.frame()
  not_permut_out <-data.frame()

  if (permut) {
    for (shuffle_index in c(1:times)) {
      # shuffle all enhancers
      all_se_temp_permut <- in_se_count_matrix
      all_se_temp_permut$S1_r1 <- sample(all_enhancer_count_matrix$S1_r1,nrow(in_se_count_matrix))
      all_se_temp_permut$S1_r2 <- sample(all_enhancer_count_matrix$S1_r2,nrow(in_se_count_matrix))
      all_se_temp_permut$S2_r1 <- sample(all_enhancer_count_matrix$S2_r1,nrow(in_se_count_matrix))
      all_se_temp_permut$S2_r2 <- sample(all_enhancer_count_matrix$S2_r2,nrow(in_se_count_matrix))

      # DESeq2 fold-change function
      fold_change <- function (permut_count_matrix) {
        permut_matrix <- permut_count_matrix[,-c(1,6:9)]
        sample_data <- data.frame(row.names = colnames(permut_matrix),
                                  condition = c("S1","S1","S2","S2"))
        dds <- DESeqDataSetFromMatrix(countData = permut_matrix,
                                      colData = sample_data,
                                      design = ~ condition)
        # calculate fold change
        dds <- DESeq(dds)
        res <- results(dds)

        # lfc shrinking
        resLFC <- lfcShrink(dds, coef="condition_S2_vs_S1", type="apeglm")
        new_resLFC <-setDT(as.data.frame(resLFC), keep.rownames = "e_merge_name")[]

        # make final output
        temp_shuffle_count_df <- merge(permut_count_matrix,new_resLFC,by="e_merge_name")

        return(temp_shuffle_count_df)
      }

      # spline fit permuted function
      fit_permut <- function(fold_change_out,merged_super_name) {
        # get new fitted value
        permut_fit <- data.frame()
        for (k in c(1:length(merged_super_name))) {
          shuffle_enhancer <- fold_change_out[which(fold_change_out$se_merge_name == merged_super_name[k]),]
          shuffle_enhancer$width_mid <- (shuffle_enhancer$start+(shuffle_enhancer$end-shuffle_enhancer$start+1)/2)/1000

          # splines bs fit
          spline_bs_fit <- lm(formula = log2FoldChange ~ bs(width_mid),data = shuffle_enhancer,
                              weights = baseMean)
          shuffle_enhancer$permut_spline_bs <- spline_bs_fit$fitted.values

          # create output dataframe with fitted value
          permut_fit <- rbind(permut_fit,shuffle_enhancer)
        }
        return(permut_fit)
      }

      # fold change of permute all enhancers in se
      all_shuffle <- fold_change(all_se_temp_permut)

      # remove 0 counts of all sample
      not_zero_shuffle <- apply(all_shuffle[,c(2:5)], 1, function(row) sum(row) != 0)
      all_shuffle <- all_shuffle[not_zero_shuffle,]

      all_permut_fit <- fit_permut(all_shuffle,merged_super_name)

      # creat output dataframe
      # all
      all_temp_new_fit_df <- as.data.frame(all_permut_fit[,c(1,16)])
      all_temp_new_fit_df$shuffle_index <- rep(shuffle_index,nrow(all_temp_new_fit_df))
      all_permut_out <- rbind(all_permut_out,all_temp_new_fit_df)

    }
  } else {
    all_permut_out <- NA
  }

  #-----------------------------------------------------------------------------
  # plot density of original bs fit and new bs fit
  if (permut) {
    plot_original <- fit_out[!is.na(fit_out$spline_bs),]
    plot_all <- all_permut_out[!is.na(all_permut_out$permut_spline_bs),]

    temp_density_plot_df <- data.frame(group = rep("original", nrow(plot_original)),
                                       w_bs_spline = plot_original$spline_bs)

    temp_density_plot_df_all <- data.frame(group = rep("All enhancer permutation", nrow(plot_all)),
                                           w_bs_spline = plot_all$permut_spline_bs)

    # all permutation in one plot
    # exclude not in SE
    density_plot_df <- rbind(temp_density_plot_df,temp_density_plot_df_all)

    # get cutoff based on slope = 1
    d <- density(plot_all$permut_spline_bs)
    slope_df <- data.frame(density = d$y[-1], log2FC = d$x[-1],
                           slope = diff(d$y)/diff(d$x))

    # remove small density values
    new_df <- slope_df[which(slope_df$density > 0.01),]
    # normalize x and y-axis
    x_unit <- max(new_df$log2FC)-min(new_df$log2FC)
    y_unit <- max(new_df$density)-min(new_df$density)
    r_xy <- y_unit/x_unit

    root_point_1 <- solveroot(slope_df$log2FC,slope_df$slope,y0 = r_xy)
    root_point_2 <- solveroot(slope_df$log2FC,slope_df$slope,y0 = -r_xy)

    fc_cutoff <- append(root_point_1,root_point_2)
    u_cut_v <- fc_cutoff[fc_cutoff > 0.5 & fc_cutoff < 2]
    l_cut_v <- fc_cutoff[fc_cutoff > -2 & fc_cutoff < -0.5]

    # check if there are number
    if (length(u_cut_v) != 0) {
      u_cut <- min(u_cut_v)
    } else if (length(l_cut_v) != 0) {
      u_cut <- -max(l_cut_v)
    } else {
      u_cut <- 1.5
    }

    if (length(l_cut_v) != 0) {
      l_cut <- max(l_cut_v)
    } else if (length(u_cut_v) != 0) {
      l_cut <- -min(u_cut_v)
    } else {
      l_cut <- -1.5
    }

    final_cutoff <- c(l_cut,u_cut)

    # plot density with cutoff line
    density_p <- ggplot(density_plot_df, aes(x = w_bs_spline,color = group))+
      geom_density()+
      geom_vline(xintercept=final_cutoff[1],linetype = "dashed",color="red")+
      geom_vline(xintercept=final_cutoff[2],linetype = "dashed",color="red")+
      theme_classic()+
      ggtitle("merged_enhancer Density")+
      xlab("Weighted lowess fitted log2FC")+
      ylab("Density")+
      xlim(c(-5,5))+
      theme(axis.text.x = element_text(angle = 0,hjust = 1),
            plot.title = element_text(hjust = 0.5))

  } else {
    density_p <- NA
    final_cutoff <- NA
  }

  # creat final output list
  out <- list()
  out$se_fit_df <- fit_out
  out$se_permutation_df <-all_permut_out
  out$cutoff <- final_cutoff
  out$density_plot <- density_p
  return(out)
}
