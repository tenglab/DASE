#' Significance estimation with permutation
#'
#' A function to estimate the significance cutoff with permutation
#'
#' @details
#' A permutation of counts of enhancers in SEs is used to decide the
#' cutoff of significant log2FC with shuffling pool. Cutoff was calculate by slope equal to range(y)/range(x).
#' The max value in range(-2,-1) and min value in range(1,2) are set to be lower and upper cutoff.
#' Coupled with smooth spline function.
#'
#' @param spline_fit_out enhancer output file from SEfitspline
#' @param sample_pool enhancer counts shuffling pool
#' (header must be "e_merge_name","S1_r1","S1_r2","S2_r1","S2_r2", "se_merge_name"))
#' @param times permutation times (default=10)
#' @param permut if you want permutation (default=TRUE)
#' @param c1_n number of replicates (samples) in condition 1.
#' @param c2_n number of replicates (samples) in condition 2.
#'
#' @return
#' se_permutation_df: permutation data frame
#' cutoff: final cutoff
#' cutoff_candidate: all cutoff candidates
#' density_plot: density plot of permuted and original log2 fold change
#'
#' @import DESeq2
#' @import Rsubread
#' @import apeglm
#'
#' @export
#' @examples
#' # default permutation 10 times
#' permut_out <- SEpermut_smooth(fit_spline_out_df,sample_pool,c1_n=2,c2_n=2)
#'
#' # permutation 5 times
#' permut_out <- SEpermut_smooth(fit_spline_out_df,sample_pool,times=5,c1_n=2,c2_n=2)
#'
#' # no permutation
#' permut_out <- SEpermut_smooth(fit_spline_out_df,sample_pool,permut=F,c1_n=2,c2_n=2)

SEpermut_smooth <- function(spline_fit_out,sample_pool,permut=T,times=10,c1_n,c2_n) {
  permut_out <- data.frame()
  if (permut) {
    for (shuffle_i in c(1:times)) {
      print(paste0("Permutation: ",shuffle_i))
      # shuffle counts
      col_list <- c('e_merge_name',colnames(spline_fit_out)[6:(5+c1_n+c2_n)],
                    'width_mid','se_merge_name')
      temp_permut <- subset(spline_fit_out,select=col_list)

      for (si in 2:(1+c1_n+c2_n)) {
        temp_permut[[si]] <- sample(sample_pool[[si]],size=nrow(temp_permut))
      }


      # DESeq2 fold-change
      # make count matrix
      count_matrix <- as.data.frame(subset(temp_permut,select=c(colnames(spline_fit_out)[6:(5+c1_n+c2_n)])))
      rownames(count_matrix) <- temp_permut$e_merge_name

      # remove row.sum = 0
      no_zero_count_matrix <- count_matrix[rowSums(count_matrix)>0,]
      sample_data <- data.frame(row.names = colnames(no_zero_count_matrix),
                                condition = c(rep("C1",c1_n),rep("C2",c2_n)))
      dds <- DESeqDataSetFromMatrix(countData = no_zero_count_matrix,
                                    colData = sample_data,
                                    design = ~ condition)

      # calculate fold change
      dds <- DESeq(dds)
      res <- results(dds)

      # save raw counts and normalized counts
      normalized_count <- counts(dds,normalized=TRUE)
      colnames(normalized_count) <- paste0(colnames(normalized_count),"_norm")
      normalized_count <- setDT(as.data.frame(normalized_count), keep.rownames = "e_merge_name")[]

      # lfc shrinking
      resLFC <- lfcShrink(dds, coef="condition_C2_vs_C1", type="apeglm")
      new_resLFC <-setDT(as.data.frame(resLFC), keep.rownames = "e_merge_name")[]

      # make final output
      temp_shuffle_count_df_temp <- merge(temp_permut,normalized_count,by="e_merge_name")
      temp_shuffle_count_df <- merge(temp_shuffle_count_df_temp,new_resLFC,by="e_merge_name")

      # add max_mean
      norm_start <-  1+c1_n+c2_n+3
      temp_shuffle_count_df$C1_mean <- rowMeans(temp_shuffle_count_df[,norm_start:(norm_start+c1_n-1)])
      temp_shuffle_count_df$C2_mean <- rowMeans(temp_shuffle_count_df[,(norm_start+c1_n):(norm_start+c1_n+c2_n-1)])
      temp_shuffle_count_df$max_mean <- apply(temp_shuffle_count_df[,c('C1_mean','C2_mean')],1,FUN=max)

      # get new fitted value
      permut_fit <- data.frame()
      se_name <- unique(temp_shuffle_count_df$se_merge_name)
      for (k in c(1:length(se_name))) {
        # print step information
        # if(k %% 100==0) {
        #   # Print on the screen some message
        #   print(paste0("SE: ",k))
        # }

        shuffle_enhancer <- temp_shuffle_count_df[which(temp_shuffle_count_df$se_merge_name == se_name[k]),]
        shuffle_enhancer$percent <- shuffle_enhancer$max_mean/sum(shuffle_enhancer$max_mean)*100
        shuffle_enhancer <- shuffle_enhancer[order(-shuffle_enhancer$percent)]
        shuffle_enhancer$cumsum <- cumsum(shuffle_enhancer$percent)

        # check if only 1 enhancer
        if (nrow(shuffle_enhancer) == 1) {
          spline_bs_fit <- lm(log2FoldChange~bs(width_mid,
                                                  degree=2),
                              data = shuffle_enhancer,
                              weights = max_mean)

          shuffle_enhancer$permut_spline_bs <- spline_bs_fit$fitted.values

          # save fitted.values for permutation
          shuffle_enhancer$permut_spline_bs <- spline_bs_fit$fitted.values

        } else {
          # get df based on number of enhancer whose cumsum is > 95%
          n_top <- min(which(shuffle_enhancer$cumsum > 95))

          # different df based on n_top
          if (n_top <= 4) {
            # df=2
            spline_bs_fit <- lm(log2FoldChange~bs(width_mid,
                                                    degree=2,
                                                    Boundary.knots = c(min(width_mid)-5,max(width_mid)+5)),
                                data = shuffle_enhancer,
                                weights = max_mean)
            # save fitted.values for permutation
            shuffle_enhancer$permut_spline_bs <- spline_bs_fit$fitted.values

          } else if (n_top >= 5) {
              # df=4
              spline_bs_fit <- smooth.spline(x=shuffle_enhancer$width_mid,y=shuffle_enhancer$log2FoldChange,
                                             w=shuffle_enhancer$max_mean,
                                             df=4)

              # save fitted.values for permutation
              shuffle_enhancer$permut_spline_bs <- spline_bs_fit$y
          }
        }

        permut_fit <- rbind(permut_fit,shuffle_enhancer)
      }

      # creat output dataframe
      temp_new_fit_df <- as.data.frame(permut_fit[,c('e_merge_name','se_merge_name',
                                                     'permut_spline_bs')])
      temp_new_fit_df$shuffle_index <- rep(shuffle_i,nrow(temp_new_fit_df))
      permut_out <- rbind(permut_out,temp_new_fit_df)
    }
  } else {
    permut_out <- NA
  }

  #-----------------------------------------------------------------------------
  # plot density of original bs fit and new bs fit
  if (permut) {
    plot_original <- spline_fit_out[!is.na(spline_fit_out$spline_bs),]
    plot_shuffle <- permut_out[!is.na(permut_out$permut_spline_bs),]

    temp_density_plot_df <- data.frame(group = rep("original", nrow(plot_original)),
                                       w_bs_spline = plot_original$spline_bs)

    temp_density_plot_df_2 <- data.frame(group = rep("permutation", nrow(plot_shuffle)),
                                         w_bs_spline = plot_shuffle$permut_spline_bs)

    density_plot_df <- rbind(temp_density_plot_df,temp_density_plot_df_2)

    # get cutoff based on slope = 1
    d <- density(plot_shuffle$permut_spline_bs)
    slope_df <- data.frame(density = d$y[-1], log2FC = d$x[-1],
                           slope = diff(d$y)/diff(d$x))

    # remove small density values
    new_df <- slope_df[which(slope_df$density > 0.02),]

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
      geom_vline(xintercept=final_cutoff[1],linetype = "dashed")+
      geom_vline(xintercept=final_cutoff[2],linetype = "dashed")+
      theme_classic()+
      ggtitle("merged_enhancer Density")+
      xlab("Weighted b-spline fitted log2FC")+
      ylab("Density")+
      xlim(c(-5,5))+

      theme(axis.text.x = element_text(angle = 0,hjust = 1),
            plot.title = element_text(hjust = 0.5))

  } else {
    density_p <- NA
    final_cutoff <- NA
    fc_cutoff <- NA
  }

  # creat final output list
  out <- list()
  out$se_permutation_df <-permut_out
  out$cutoff <- final_cutoff
  out$cutoff_candidate <- fc_cutoff
  out$density_plot <- density_p
  return(out)
}

