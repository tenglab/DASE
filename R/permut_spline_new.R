# shuffle counts of each enhancer n times (default 10), got new fold change,

permut_spline_new <- function(spline_fit_out,sample_pool,permut=T,times=10) {

  permut_out <- data.frame()
  if (permut) {
    for (shuffle_i in c(1:times)) {
      print(paste0("Permutation: ",shuffle_i))
      # shuffle counts
      temp_permut <- spline_fit_out[,c('e_merge_name','S1_r1','S1_r2','S2_r1','S2_r2',
                                       'width_mid','baseMean.1','se_merge_name')]
      temp_permut$S1_r1 <- sample(sample_pool$S1_r1,size=nrow(temp_permut))
      temp_permut$S1_r2 <- sample(sample_pool$S1_r2,size=nrow(temp_permut))
      temp_permut$S2_r1 <- sample(sample_pool$S2_r1,size=nrow(temp_permut))
      temp_permut$S2_r2 <- sample(sample_pool$S2_r2,size=nrow(temp_permut))

      # DESeq2 fold-change
      # make count matrix
      count_matrix <- as.data.frame(temp_permut[,c('S1_r1','S1_r2','S2_r1','S2_r2')])
      rownames(count_matrix) <- temp_permut$e_merge_name

      # remove row.sum = 0
      no_zero_count_matrix <- count_matrix[rowSums(count_matrix)>0,]
      sample_data <- data.frame(row.names = colnames(no_zero_count_matrix),
                                condition = c("S1","S1","S2","S2"))
      dds <- DESeqDataSetFromMatrix(countData = no_zero_count_matrix,
                                    colData = sample_data,
                                    design = ~ condition)

      # calculate fold change
      dds <- DESeq(dds)
      res <- results(dds)

      # save raw counts and normalized counts
      normalized_count <- counts(dds,normalized=TRUE)
      colnames(normalized_count) <- c("S1_r1_norm","S1_r2_norm","S2_r1_norm","S2_r2_norm")
      normalized_count <- setDT(as.data.frame(normalized_count), keep.rownames = "e_merge_name")[]

      # lfc shrinking
      resLFC <- lfcShrink(dds, coef="condition_S2_vs_S1", type="apeglm")
      new_resLFC <-setDT(as.data.frame(resLFC), keep.rownames = "e_merge_name")[]

      # make final output

      temp_shuffle_count_df_temp <- merge(temp_permut,normalized_count,by="e_merge_name")
      temp_shuffle_count_df <- merge(temp_shuffle_count_df_temp,new_resLFC,by="e_merge_name")

      # add max_mean
      temp_shuffle_count_df$s1_mean <- (temp_shuffle_count_df$S1_r1_norm+temp_shuffle_count_df$S1_r2_norm)/2
      temp_shuffle_count_df$s2_mean <- (temp_shuffle_count_df$S2_r1_norm+temp_shuffle_count_df$S2_r2_norm)/2
      temp_shuffle_count_df$max_mean <- apply(temp_shuffle_count_df[,c('s1_mean','s2_mean')],1,FUN=max)

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

          } else if (n_top >= 6) {
              # df=4
              spline_bs_fit <- lm(log2FoldChange~bs(width_mid,
                                                      degree=4,
                                                      Boundary.knots = c(min(width_mid)-5,max(width_mid)+5)),
                                  data = shuffle_enhancer,
                                  weights = max_mean)

          } else {
            # df=3
            spline_bs_fit <- lm(log2FoldChange~bs(width_mid,
                                                    degree=3,
                                                    Boundary.knots = c(min(width_mid)-5,max(width_mid)+5)),
                                data = shuffle_enhancer,
                                weights = max_mean)

          }

          # save fitted.values for permutation
          shuffle_enhancer$permut_spline_bs <- spline_bs_fit$fitted.values

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
  }

  # creat final output list
  out <- list()
  out$se_permutation_df <-permut_out
  out$cutoff <- final_cutoff
  out$cutoff_candidate <- fc_cutoff
  out$density_plot <- density_p
  return(out)
}

