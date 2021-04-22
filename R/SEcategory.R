#' Category and ranking of SE
#'
#' A function to get the category of SEs based on SE segment percentage file
#'
#' @details
#' This function will create the category of SEs (similar, strengthen/weaken, shorten,
#' shifting, V-shape, other) based on the segment percentage. Within each category, SE are
#' firstly ranked by segment percentage, then by average fold change of enhancers within
#'
#' @param se_seg_df SE segment percentage file form SEpattern
#' @param e_deseq_out enhancer deseq out file from enhancerFoldchange
#'
#' @return
#' A dataset of final SE category and ranking.
#'
#' @export
#' @examples
#' catRank_list <- SEcategory(segment_percentage_df,enhancerFoldchange_out)
#'

SEcategory <- function(se_seg_df,e_deseq_out) {

  # remove all segement whose percentage is less than 0.05
  se_seg_filter <- se_seg_df[which(se_seg_df$seg_percent > 0.05),]
  se_name <- unique(se_seg_df$se_merge_name)

  #---------------------------------------------------------------------------
  # get SE category
  #---------------------------------------------------------------------------
  se_seg_w_cat_df <- data.frame()
  for (se_index in c(1:length(se_name))) {
    temp_seg_df <- se_seg_filter[which(se_seg_filter$se_merge_name == se_name[se_index]),]
    temp_seg_df <- temp_seg_df[order(temp_seg_df$seg_name),]
    n_row <- nrow(temp_seg_df)

    # calculate new percentage after remove the small one
    total_new <- sum(temp_seg_df$seg_percent)
    temp_seg_df$seg_percent <- temp_seg_df$seg_percent/total_new

    #---------------------------------------------------------------------------
    # first merge the segment with same location in a row for (nrow > 1)
    #---------------------------------------------------------------------------
    if (n_row > 1) {
      # check if segment location is same in a row
      flag_list <- c()
      for (seg_index in c(1:(n_row-1))) {
        if (temp_seg_df$seg_loc[seg_index] == temp_seg_df$seg_loc[seg_index+1]) {
          flag <- T
          flag_list <- append(flag_list,flag)
        } else {
          flag <- F
          flag_list <- append(flag_list,flag)
        }
      }

      # if there are same location in a row, merge them
      if (any(flag_list == T)) {
        # find which segment are the TRUE
        merge_start_point <- which(flag_list == T)

        # change seg_name based on where TRUE values are
        if (length(merge_start_point) > 1) {
          merge_start_point <- append(merge_start_point,max(merge_start_point)+1)
          for (index in c(1:(length(merge_start_point)-1))) {
            if (merge_start_point[index] == (merge_start_point[index+1]-1)) {
              temp_seg_df$seg_name[merge_start_point[index]] <- temp_seg_df$seg_name[max(merge_start_point)]
            } else {
              temp_seg_df$seg_name[merge_start_point[index]] <- temp_seg_df$seg_name[merge_start_point[index]+1]
            }
          }
        } else {
          temp_seg_df$seg_name[merge_start_point] <- temp_seg_df$seg_name[merge_start_point+1]
        }

        # merge common segment
        temp_cat_df <- aggregate(data=temp_seg_df,seg_percent~.,sum)
        temp_cat_df <- temp_cat_df[order(temp_cat_df$seg_name),]
      } else {
        temp_cat_df <- temp_seg_df
      }
    } else {
      temp_cat_df <- temp_seg_df
    }

    # add new number of segments
    temp_cat_df$seg_number <- rep(nrow(temp_cat_df),nrow(temp_cat_df))
    #---------------------------------------------------------------------------
    # then get catogory
    #---------------------------------------------------------------------------

    # initial variables
    temp_cat_df$category <- NA
    temp_cat_df$direction <- NA

    # check only 1 segment (middle, upper, or lower)
    if (unique(temp_cat_df$seg_number) == 1) {
      if (temp_cat_df$seg_loc == "mid") {
        temp_cat_df$category <- rep("Similar",nrow(temp_cat_df))
        temp_cat_df$direction <- "none"

      } else {
        temp_cat_df$category <- rep("Strengthen/weaken",nrow(temp_cat_df))
        if (temp_cat_df$seg_loc == "lower") {
          temp_cat_df$direction <- "-"
        } else {
          temp_cat_df$direction <- "+"
        }
      }
    }

    #------------------------------------------------------------------------------
    # check 2 segment
    if (unique(temp_cat_df$seg_number) == 2) {

      # extract mid and other segments
      temp_mid <- temp_cat_df[which(temp_cat_df$seg_loc == "mid"),]
      temp_other <- temp_cat_df[which(temp_cat_df$seg_loc != "mid"),]

      if (temp_mid$seg_percent >= 0.9) {
        temp_cat_df$category <- rep("Similar",nrow(temp_cat_df))
        temp_cat_df$direction <- "none"

      } else if (temp_mid$seg_percent < 0.9 & temp_mid$seg_percent >= 0.5) {
        temp_cat_df$category <- rep("Shorten",nrow(temp_cat_df))
        if (temp_other$seg_loc == "lower") {
          temp_cat_df$direction <- "-"
        } else {
          temp_cat_df$direction <- "+"
        }

      } else {
        temp_cat_df$category <- rep("Strengthen/weaken",nrow(temp_cat_df))
        if (temp_other$seg_loc == "lower") {
          temp_cat_df$direction <- "-"
        } else {
          temp_cat_df$direction <- "+"
        }
      }
    }

    #------------------------------------------------------------------------------
    # check 3 segment
    if (unique(temp_cat_df$seg_number) == 3) {

      # extract each segments
      temp_p_1 <- temp_cat_df[1,]
      temp_p_2 <- temp_cat_df[2,]
      temp_p_3 <- temp_cat_df[3,]

      # check shifting
      if (length(unique(temp_cat_df$seg_loc)) == 3) {
        temp_cat_df$category <- rep("Shifting",nrow(temp_cat_df))
        if (temp_p_1$seg_loc == "lower") {
          temp_cat_df$direction <- "r"
        } else {
          temp_cat_df$direction <- "l"
        }

      # check only have upper/mid or lower/mid and part_2 is mid
      } else if ( length(unique(temp_cat_df$seg_loc)) == 2
                  & temp_p_2$seg_loc == "mid") {
          if ( temp_p_2$seg_percent >= 0.9) {
            temp_cat_df$category <- rep("Similar",nrow(temp_cat_df))
          } else if (temp_p_2$seg_percent < 0.9 & temp_p_2$seg_percent >= 0.5) {
              temp_cat_df$category <- rep("Shorten",nrow(temp_cat_df))
              if (temp_p_1$seg_loc == "lower") {
                temp_cat_df$direction <- "-"
              } else {
                temp_cat_df$direction <- "+"
              }
          } else {
            temp_cat_df$category <- rep("Strengthen/weaken",nrow(temp_cat_df))
            if (temp_p_1$seg_loc == "lower") {
              temp_cat_df$direction <- "-"
            } else {
              temp_cat_df$direction <- "+"
            }
          }

      # check part_2 is lower or upper
      } else if ( length(unique(temp_cat_df$seg_loc)) == 2 &
                  (temp_p_2$seg_loc == "lower" |
                   temp_p_2$seg_loc == "upper") ) {
        if ( temp_p_2$seg_percent >= 0.5) {
          temp_cat_df$category <- rep("Strengthen/weaken",nrow(temp_cat_df))
          if (temp_p_2$seg_loc == "lower") {
            temp_cat_df$direction <- "-"
          } else {
            temp_cat_df$direction <- "+"
          }
        } else if (temp_p_2$seg_percent < 0.1) {

          temp_cat_df$category <- rep("Similar",nrow(temp_cat_df))
          temp_cat_df$direction <- "none"

        } else if (temp_p_2$seg_percent < 0.5 &
                   temp_p_2$seg_percent > 0.1 &
                   abs(temp_p_1$seg_percent - temp_p_3$seg_percent) < 0.2) {
          temp_cat_df$category <- rep("V_shape",nrow(temp_cat_df))
          if (temp_p_2$seg_loc == "lower") {
            temp_cat_df$direction <- "-"
          } else {
            temp_cat_df$direction <- "+"
          }

        } else if (temp_p_2$seg_percent < 0.5 &
                   temp_p_2$seg_percent > 0.1 &
                   abs(temp_p_1$seg_percent - temp_p_3$seg_percent) > 0.2) {
          temp_cat_df$category <- rep("Shorten",nrow(temp_cat_df))
          if (temp_p_2$seg_loc == "lower") {
            temp_cat_df$direction <- "-"
          } else {
            temp_cat_df$direction <- "+"
          }

        } else {
          temp_cat_df$category <- rep("Other",nrow(temp_cat_df))
          temp_cat_df$direction <- "none"
        }
      }
    }

    #------------------------------------------------------------------------------
    # more than 4 segments
    if (unique(temp_cat_df$seg_number) >= 4) {
      temp_cat_df$category <- rep("Other",nrow(temp_cat_df))
      temp_cat_df$direction <- "none"
    }

    #------------------------------------------------------------------------------

    se_seg_w_cat_df <- rbind(se_seg_w_cat_df,temp_cat_df)
  }

  #---------------------------------------------------------------------------
  # rank SE within each category
  #---------------------------------------------------------------------------
  # get average log2Foldchange for each SE

  mean_fc_df <- data.frame()
  for (se in c(1:length(se_name))) {
    temp_se <- e_deseq_out[which(e_deseq_out$se_merge_name == se_name[se]),]
    temp_se$mean_FC <- rep(mean(temp_se$log2FoldChange.1),nrow(temp_se))
    temp_out <- unique(temp_se[,c(25,26)])
    mean_fc_df <- rbind(mean_fc_df,temp_out)
  }

  # rank similar
  temp_similar_df <- se_seg_w_cat_df[which(se_seg_w_cat_df$category=="Similar" & se_seg_w_cat_df$seg_loc == "mid"),]
  similar_df <- aggregate(data=temp_similar_df[,-5],seg_percent~.,sum)

  similar_merge_df <- merge(similar_df,mean_fc_df,by="se_merge_name")
  rank_sim_df <- similar_merge_df[order(-similar_merge_df$seg_percent,abs(similar_merge_df$mean_FC)),]
  rank_sim_df$rank <- seq(1:nrow(rank_sim_df))

  # rank strength/weaken
  temp_strength_df <- se_seg_w_cat_df[which(se_seg_w_cat_df$category=="Strengthen/weaken"& se_seg_w_cat_df$seg_loc != "mid"),]
  strength_df <- aggregate(data=temp_strength_df[,-5],seg_percent~.,sum)

  strength_merge_df <- merge(strength_df,mean_fc_df,by="se_merge_name")
  rank_str_df <- strength_merge_df[order(-strength_merge_df$seg_percent,-abs(strength_merge_df$mean_FC)),]
  rank_str_df$rank <- seq(1:nrow(rank_str_df))

  # rank shorten
  temp_short_df <- se_seg_w_cat_df[which(se_seg_w_cat_df$category=="Shorten"& se_seg_w_cat_df$seg_loc == "mid"),]
  short_df <- aggregate(data=temp_short_df[,-5],seg_percent~.,sum)

  short_merge_df <- merge(short_df,mean_fc_df,by="se_merge_name")
  rank_short_df <- short_merge_df[order(short_merge_df$seg_percent,-abs(short_merge_df$mean_FC)),]
  rank_short_df$rank <- seq(1:nrow(rank_short_df))

  # rank V-shape
  temp_v_df <- se_seg_w_cat_df[which(se_seg_w_cat_df$category=="V_shape"& se_seg_w_cat_df$seg_loc != "mid"),]
  v_df <- aggregate(data=temp_v_df[,-5],seg_percent~.,sum)

  v_merge_df <- merge(v_df,mean_fc_df,by="se_merge_name")
  rank_v_df <- v_merge_df[order(-v_merge_df$seg_percent,-abs(v_merge_df$mean_FC)),]
  rank_v_df$rank <- seq(1:nrow(rank_v_df))

  # rank Shifting
  temp_shift_df <- se_seg_w_cat_df[which(se_seg_w_cat_df$category=="Shifting"& se_seg_w_cat_df$seg_loc == "mid"),]

  shift_merge_df <- merge(temp_shift_df,mean_fc_df,by="se_merge_name")
  rank_shift_df <- shift_merge_df[order(shift_merge_df$seg_percent,-abs(shift_merge_df$mean_FC)),]
  rank_shift_df$rank <- seq(1:nrow(rank_shift_df))

  # other no rank
  other_shift_df <- se_seg_w_cat_df[which(se_seg_w_cat_df$category=="Other"),]

  rank_other_df <- unique(other_shift_df[,-c(5:7)])
  rank_other_df$rank <- rep("none",nrow(rank_other_df))

  # make output file
  se_cat_rank <- unique(rbind(rank_sim_df[,-c(5,8,9)],rank_str_df[,-c(5,8,9)],
                       rank_sim_df[,-c(5,8,9)],rank_short_df[,-c(5,8,9)],
                       rank_v_df[,-c(5,8,9)],rank_shift_df[,-c(5:7,10)],rank_other_df))
  return(se_cat_rank)
}


