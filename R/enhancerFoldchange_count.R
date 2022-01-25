#' Calculate the fold change of enhancers
#'
#' A function to calculate the fold change of enhancers base on two conditions.
#'
#' @details
#' This function will extract the enhancers which in the range of SEs with count table of pooled enhancers.
#' Then enhancers are merged with no gaps.
#' DESeq2 is used to calculate the fold change of those enhancers
#'
#' @param e_df pooled macs2 output file *_peaks.narrowPeak of all samples.
#' @param se_df merged SE metadata from SEfilter.
#' @param raw_count raw count matrix of pooled enhancers in each sample.
#'                  The order of sample columns need to be "condition1_rep1, condition1_rep2,...,condition2_rep1,condition2_rep2"
#' @param c1_n number of replicates (samples) in condition 1.
#' @param c2_n number of replicates (samples) in condition 2.
#'
#' @return
#' A list of 2 datasets: enhancer foldchang results dataset and sizefactor used in DESeq2
#'
#' @import data.table
#' @import DESeq2
#' @import Rsubread
#' @import apeglm
#' @import tidyr
#'
#' @export
#' @examples
#' foldchange_list <- enhancerFoldchange_count(e_df,se_df,raw_count,c1_n=2,c2_n=2)
#'

enhancerFoldchange_count <- function(e_df,se_df,raw_count,c1_n,c2_n) {
  colnames(raw_count)[1] <- "enhancer"

  # extract chr, start, end from raw_count
  raw_count_temp<- separate(data = raw_count, col = enhancer, into = c("chr", "start","end"), sep = "_")

  # merge enhancers
  ir <- IRanges(as.integer(raw_count_temp$start),
                as.integer(raw_count_temp$end),
                names = raw_count_temp$chr)

  e_merge_by_chr <- rbindlist(lapply(split(ir, names(ir)),
                                     function(x) as.data.table(reduce(x,min.gapwidth = 0))),
                              idcol = "chr")
  colnames(e_merge_by_chr) <- c("chr","start","end","width")

  # create merged enhancer names
  e_merge_by_chr$e_merge_name <- paste0(e_merge_by_chr$chr,"_",e_merge_by_chr$start,"_",e_merge_by_chr$end)

  # find overlaps with new merged enhancers and aggregate the count of them
  ir_merge <- IRanges(e_merge_by_chr$start,
                      e_merge_by_chr$end,
                      names=e_merge_by_chr$chr)

  hit <- findOverlaps(ir,ir_merge)

  raw_count_temp$e_merge_name[queryHits(hit)] <- e_merge_by_chr$e_merge_name[subjectHits(hit)]

  # aggregate count
  count_matrix <- aggregate(.~e_merge_name,data=raw_count_temp[-c(1:3)],sum)


  # DESeq2 fold-change for enhancer in SE only
  # make count matrix
  row.names(count_matrix) <- count_matrix$e_merge_name
  count_matrix <- count_matrix[,-1]

  # remove enhancers with 0 value for all sample
  not_zero_in_se <- apply(count_matrix, 1, function(row) sum(row) != 0)
  count_matrix <- count_matrix[not_zero_in_se,]

  # change column names
  colnames(count_matrix) <- gsub("S","C",colnames(count_matrix))
  # make sample data for DESeq2
  sample_data <- data.frame(row.names = colnames(count_matrix),
                            condition = c(rep("C1",c1_n),rep("C2",c2_n)))

  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = sample_data,
                                design = ~ condition)

  # run DEseq2
  dds <- DESeq(dds)
  res <- results(dds)

  # save raw counts and normalized counts

  raw_count <- counts(dds)
  normalized_count <- counts(dds,normalized=TRUE)
  colnames(normalized_count) <- paste0(colnames(normalized_count),"_norm")

  # shrinking LFC and MA plot
  resLFC <- lfcShrink(dds, coef="condition_C2_vs_C1", type="apeglm")

  # make final output
  temp_output <- data.frame()
  temp_output <- as.data.frame(cbind(raw_count,normalized_count,res,resLFC))
  temp_output <- setDT(temp_output, keep.rownames = "e_merge_name")[]

  # combine DESeq2 results with merged enhancer profile
  merged_out_df <- merge(e_merge_by_chr,temp_output,by="e_merge_name")

  # assign merged super-enhancer name to merged enhancer deseq outfile
  final_out_df <- data.frame()
  for (i in c(1:nrow(se_df))) {
    temp_merge <- merged_out_df[which(merged_out_df$chr==se_df$chr[i] &
                                        merged_out_df$start>=se_df$start[i] &
                                        merged_out_df$end<= se_df$end[i]),]
    temp_merge$se_merge_name <- rep(se_df$se_merge_name[i],nrow(temp_merge))
    final_out_df <- rbind(final_out_df,temp_merge)
  }

  # save size_factor
  size_factor <- sizeFactors(dds)

  not_in_se_deseq <- temp_output[!(temp_output$e_merge_name %in% final_out_df$e_merge_name),]

  # make output list
  out <- list()
  out$lfc_shrink <- resLFC
  out$enhancer_deseq_result <- final_out_df
  out$size_factor <- size_factor
  out$not_in_se_deseq <- not_in_se_deseq
  return(out)
}
