#' Filter super-enhancer
#'
#' A function to filter super-enhancers (SEs)
#'
#' @details
#' This function firstly remove SEs whose width are within lower and upper 1 percent.
#' Then filter out SEs in the range of super-enhancer blacklist.
#' Finally, remained SEs from different samples are merged to the longest representative SE.
#'
#' @param se_in merged SE region bed file of all samples.
#' @param bl_file super-enhancer blacklist bed file download from ENCODE (ENCFF356LFX).
#' @param custom_range a vector of extra customized blacklist to ignore.
#' Format: c(chr:start-end, chr:start-end, ...).
#'
#' @return
#' se_filtered_no_bl: filtered SE datasets with merged SE names,
#' se_filtered_in_bl: SE datasets identified in the blacklist
#' se_merged_meta: merged SE metadata
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

SEfilter <- function (se_in,bl_file,custom_range) {
  # get width of each SE
  se_in$width <- se_in$V3-se_in$V2+1

  # remove SEs whose width are within lower and upper 1%
  se_sort_df <- se_in[order(se_in$width),]
  se_filter_df <- se_sort_df[c((ceiling(nrow(se_in)*0.01)+1):
                                (nrow(se_in)-ceiling(nrow(se_in)*0.01)+1)),]

  # check if a SE is in blacklist
  if (missing(custom_range) & !missing(bl_file)) {
    new_bl_file <- bl_file
    new_bl_file$V2 <- as.integer(new_bl_file$V2)
    new_bl_file$V3 <- as.integer(new_bl_file$V3)
  } else if (!missing(custom_range) & !missing(bl_file)) {
    extra_bl <- matrix(unlist(strsplit(custom_range,":|-")),
                       ncol=3, byrow=TRUE)
    new_bl_file <- rbind(bl_file,extra_bl)
    new_bl_file$V2 <- as.integer(new_bl_file$V2)
    new_bl_file$V3 <- as.integer(new_bl_file$V3)
  } else if (!missing(custom_range) & missing(bl_file)){
    extra_bl <- matrix(unlist(strsplit(custom_range,":|-")),
                       ncol=3, byrow=TRUE)
    new_bl_file <- data.frame(V1=extra_bl[,1],
                              V2=extra_bl[,2],
                              V3=extra_bl[,3])
    new_bl_file$V2 <- as.integer(new_bl_file$V2)
    new_bl_file$V3 <- as.integer(new_bl_file$V3)

  } else if (missing(custom_range) & missing(bl_file)){
    new_bl_file <- data.frame()
  }

  if (!missing(bl_file) | !missing(custom_range)) {
    se_filter_df_in_bl <-data.frame()
    for (i in c(1:nrow(new_bl_file))) {
      in_bl <- se_filter_df[which(se_filter_df$V1 == new_bl_file$V1[i] &
                                    se_filter_df$V2 >= new_bl_file$V2[i] &
                                    se_filter_df$V3  <= new_bl_file$V3[i]),]
      se_filter_df_in_bl <- rbind(se_filter_df_in_bl,in_bl)
      se_filter_df_no_bl <- se_filter_df[!(se_filter_df$V4 %in% se_filter_df_in_bl$V4),]

    }
  } else {
    se_filter_df_no_bl <- se_filter_df
    se_filter_df_in_bl <- NA
  }

  # merge SEs from each sample to the longest representative SE
  ir <- IRanges(se_filter_df_no_bl$V2,
                se_filter_df_no_bl$V3,
                names=se_filter_df_no_bl$V1)

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
