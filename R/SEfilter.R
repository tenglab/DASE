#' Filter super-enhancer
#'
#' A function to filter super-enhancers (SEs)
#'
#' @details
#' This function firstly remove SEs whose width are within lower and upper 1 percent.
#' Then filter out SEs in the range of super-enhancer blacklist.
#' Finally, remained SEs from different samples are merged to the longest representative SE.
#'
#' @param se_in pooled ROSE output file *_peaks_Gateway_SuperEnhancers.bed of all samples.
#' Header needs to be "CHROM,START,STOP,REGION_ID,Signal,STRAND"
#' @param bl_file super-enhancer blacklist bed file download from ENCODE (ENCFF356LFX) (default=FALSE)
#' @param has_bl_file if there is a blacklist file (default=FALSE)
#'
#' @return
#' A list of 3 datasets: filtered SE datasets with merged SE names,
#' SE datasets identified in the blacklist, and merged SE metadata
#'
#' @import data.table
#'
#' @export
#' @examples
#' # run without blacklist file
#' se_list <- SEfilter(se_in)
#'
#' # run with blacklist file
#' se_list <- SEfilter(se_in,bl_file,has_bl_file=TRUE)
#'

SEfilter <- function (se_in,bl_file=FALSE,has_bl_file=FALSE) {
  # get width of each SE
  se_in$width <- se_in$STOP-se_in$START+1

  # remove SEs whose width are within lower and upper 1%
  se_sort_df <- se_in[order(se_in$width),]
  se_filter_df <- se_sort_df[c((ceiling(nrow(se_in)*0.01)+1):
                                (nrow(se_in)-ceiling(nrow(se_in)*0.01)+1)),]

  # check if a SE is in blacklist

  if (has_bl_file == TRUE) {
    se_filter_df_in_bl <-data.frame()
    for (i in c(1:nrow(bl_file))) {
      in_bl <- se_filter_df[which(se_filter_df$CHROM == bl_file$V1[i] &
                                    se_filter_df$START >= bl_file$V2[i] &
                                    se_filter_df$STOP  <= bl_file$V3[i]),]
      se_filter_df_in_bl <- rbind(se_filter_df_in_bl,in_bl)
      se_filter_df_no_bl <- se_filter_df[!(se_filter_df$REGION_ID %in% se_filter_df_in_bl$REGION_ID),]

    }
  } else {
    se_filter_df_no_bl <- se_filter_df
    se_filter_df_in_bl <- NA
  }

  # remove SE which are in the blacklist


  # merge SEs from each sample to the longest representative SE
  ir <- IRanges(se_filter_df_no_bl$START,
                se_filter_df_no_bl$STOP,
                names=se_filter_df_no_bl$CHROM)

  se_merge_by_chr <- rbindlist(lapply(split(ir, names(ir)),
                                      function(x) as.data.table(reduce(x))),
                                idcol = "chr")

  # create merged SE names
  se_merge_by_chr$se_merge_name <- apply(se_merge_by_chr[,c(1:3)],1, paste,collapse ="_" )
  se_merge_by_chr$se_merge_name <- gsub(" ","",se_merge_by_chr$se_merge_name)


  # make output list
  out <- list()
  out$se_filtered_no_bl <- se_filter_df_no_bl
  out$se_filtered_in_bl <- se_filter_df_in_bl
  out$se_merged_meta <- se_merge_by_chr
  return(out)
}
