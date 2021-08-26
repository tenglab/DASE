#' Pattern of SE
#'
#' A function to get the pattern of SEs based on enhancer fitted line
#'
#' @details
#' This function will generate plot patterns of each SE and based on the intersection of
#' cutoff and fitted line, patterns are seperated into several segments. The percentage of
#' each segments are calculated by basemean of enhancer/sum of basemean within each SE
#'
#' @param se_fit_df fitted output file from SEfitspline
#' @param cutoff_vector cutoff vector including lower and upper cutoff from SEfitspline
#' (default is c(-1.5,1.5))
#'
#' @return
#' A list of 2 datasets: plot list contain plots' object of each SE and
#' SE segment percentage dataset
#'
#' @import ggplot2
#'
#' @export
#' @examples
#' seg_list <- SEpattern(SEfitspline_out,SEfitspline_cutoff)
#'

SEpattern_demo <- function(se_fit_df,pattern_plot_df,cutoff_vector) {

  # inital output list
  plot_list <- list()
  se_segment <- data.frame()
  out <- list()

  # set significant fold change cutoffs
  l_cut <- min(cutoff_vector)
  u_cut <- max(cutoff_vector)

  # set cutoff brand to filter cutoff boundary points
  l_band <- c(l_cut*1,l_cut*1)
  u_band <- c(u_cut*1,u_cut*1)
  # loop by super-enhancer names
  se_name <- unique(se_fit_df$se_merge_name)
  for (se_index in c(1:length(se_name))) {

    # print step information
    if(se_index %% 500==0) {
      # Print on the screen some message
      print(paste0("SE: ",se_index))
    }

    #--------------------------------------------------------------
    # plot fit patterns for each SE
    #--------------------------------------------------------------
    # make plot dataframe
    temp_df <- se_fit_df[which(se_fit_df$se_merge_name==se_name[se_index]),]
    temp_plot_df <- pattern_plot_df[which(pattern_plot_df$se_merge_name==se_name[se_index]),]

    plot_title_name <- temp_plot_df$se_merge_name

    top_points <- temp_df[which(temp_df$cumsum <= 96),]
    # plot
    pattern_plot <- ggplot(temp_plot_df, aes(x=new_x, y=pred_y)) +
                          geom_point(data=temp_df, aes(x=width_mid,y=log2FoldChange.1))+
                          geom_point(data=top_points, aes(x=width_mid,y=log2FoldChange.1),col="red")+
                          geom_line()+
                          geom_hline(yintercept=l_cut,linetype = "dashed")+
                          geom_hline(yintercept=u_cut,linetype = "dashed")+
                          theme_classic()+
                          ggtitle(plot_title_name)+
                          xlab("Enhancer middle location")+
                          ylab("Log2FC")+
                          theme(plot.title = element_text(hjust = 0.5),
                                axis.text.x = element_text(angle = 45,hjust = 1),
                                axis.line.x = element_line())
    # append each plot to output plots list
    plot_list[[se_index]]<- pattern_plot

    #--------------------------------------------------------------
    # calculate porpotion of segments based on max_mean
    #--------------------------------------------------------------
    # get intersection points
    temp_x <- temp_df$width_mid
    temp_y <- temp_df$spline_bs
    x_cut_upper <- solveroot(temp_df$width_mid,temp_df$spline_bs,y0=u_band[1])
    x_cut_lower <- solveroot(temp_df$width_mid,temp_df$spline_bs,y0=l_band[1])
    x_intersection <- sort(append(x_cut_lower,x_cut_upper))
    min_width_mid <- min(temp_df$width_mid)
    max_width_mid <- max(temp_df$width_mid)
    x_cut_all <- c(min_width_mid,x_intersection,max_width_mid)

    # total
    total_width <- sum(temp_df$max_mean)

    n <- length(x_cut_all)

    #create segment name
    part_name_l <- c()

    for (p_i in c(1:(n-1))) {
      p_name <- paste("seg_",p_i,sep = "")
      part_name_l <- append(part_name_l,p_name)
    }

    # check number of width cutoff number
    part_res <- c()
    part_loc_res <- c()

    # check n == 2, means no interaction of fitted line and cutoff lines
    #
    if (n == 2) {
      max_basemean <- temp_df[which.max(temp_df$max_mean),]
      if (max_basemean$spline_bs < l_band[1]) {
        temp_lower_width <- total_width/total_width
        part_res <- append(part_res,temp_lower_width)
        part_loc_res <- append(part_loc_res,"lower")
      } else if( max_basemean$spline_bs >= l_band[1] &
                 max_basemean$spline_bs <= u_band[1] ) {
        temp_mid_width <- total_width/total_width
        part_res <- append(part_res,temp_mid_width)
        part_loc_res <- append(part_loc_res,"mid")
      } else if ( max_basemean$spline_bs >= u_band[1]) {
        temp_upper_width <- total_width/total_width
        part_res <- append(part_res,temp_upper_width,)
        part_loc_res <- append(part_loc_res,"upper")
      }
    } else {
      # have interaction of fitted line and cutoff lines
      for (m_i in c(1:(n-1))) {
        # get data point within each segment
        temp_seg <- temp_df[which(temp_df$width_mid >= x_cut_all[m_i] & temp_df$width_mid <= x_cut_all[m_i+1]),]

        if (nrow(temp_seg) != 0) {
          # get segment percentage
          seg_percent <- sum(temp_seg$max_mean)/total_width

          # using LFC of data point which has largest basemean to define segment location
          max_basemean <- temp_seg[which.max(temp_seg$max_mean),]
          if (max_basemean$spline_bs < l_band[1]) {
            part_res <- append(part_res,seg_percent,)
            part_loc_res <- append(part_loc_res,"lower")
          } else if (max_basemean$spline_bs >= l_band[1] & max_basemean$spline_bs <= u_band[1]) {
            part_res <- append(part_res,seg_percent,)
            part_loc_res <- append(part_loc_res,"mid")
          } else {
            part_res <- append(part_res,seg_percent,)
            part_loc_res <- append(part_loc_res,"upper")
          }
        } else {
          # chose cloest point as representive
          seg_percent <- 0
          part_res <- append(part_res,seg_percent,)
          part_loc_res <- append(part_loc_res,"mid")
        }
      }
    }

    # creat out pattern dataframe
    temp_per_df <- data.frame( se_merge_name = rep(se_name[se_index],length(part_name_l)),
                               total_width = total_width,
                               number_enhancer = nrow(temp_df),
                               seg_number = n-1,
                               seg_name = part_name_l,
                               seg_loc = part_loc_res,
                               seg_percent = round(part_res,digit=3)
                             )

    se_segment <- rbind(se_segment,temp_per_df)
  }

  # make final output
  out$plots <- plot_list
  out$se_segment_percent <- se_segment
  return(out)

}
