#' Calculate the fold change of enhancers
#'
#' A function to calculate the fold change of enhancers base on two conditions.
#'
#' @details
#' This function will extract the enhancers which in the range of SEs.
#' Then enhancers are merged if there gaps are within 500bps.
#' Feature count is used if provide bam file
#' rtracklayer is used if privide bw file
#' DESeq2 are used to calculate the fold change of those enhancers
#'
#' @param e_df pooled macs2 output file *_peaks.narrowPeak of all samples.
#' Header needs to be "chr, start,end,name,score,strand,signalValue,pValue,qValue,peak"
#'
#' @param se_df merged SE metadata from SEfilter
#' @param s1_r1_bw path of sample 1 replicate 1 bw file
#' @param s1_r2_bw path of sample 1 replicate 2 bw file
#' @param s2_r1_bw path of sample 2 replicate 1 bw file
#' @param s2_r2_bw path of sample 2 replicate 2 bw file
#'
#' @return
#' A list of 2 datasets: enhancer foldchang results dataset and sizefactor used in DESeq2
#'
#' @import rtracklayer
#' @import data.table
#' @import DESeq2
#' @import Rsubread
#' @import apeglm
#'
#' @export
#' @examples
#' foldchange_list <- enhancerFoldchange(pool_enhancer_df,se_meta,s1_r1_bw,s1_r2_bw,s2_r1_bw,s2_r2_bw)
#'

enhancerFoldchange_bw <- function(e_df,se_df,
                                  s1_r1_bw,s1_r2_bw,s2_r1_bw,s2_r2_bw) {

  # extract enhancers within merged super-enhancer
  in_se_df <- data.frame()
  for (i in c(1:nrow(se_df))) {

    in_se_temp <- e_df[which(e_df$chr==se_df$chr[i] &
                               e_df$start>=se_df$start[i] &
                               e_df$end<= se_df$end[i]),]
    in_se_df <- rbind(in_se_df,in_se_temp)
  }

  # merge enhancers with gaps less than 500bps
  ir <- IRanges(in_se_df$start,
                in_se_df$end,
                names = in_se_df$chr)

  e_merge_by_chr <- rbindlist(lapply(split(ir, names(ir)),
                                     function(x) as.data.table(reduce(x,min.gapwidth = 500))),
                              idcol = "chr")

  # create merged enhancer names
  e_merge_by_chr$e_merge_name <- apply(e_merge_by_chr[,c(1:3)],1, paste,collapse ="_" )
  e_merge_by_chr$e_merge_name <- gsub(" ","",e_merge_by_chr$e_merge_name)

  # get counts from bw file for one sample

  # function to count bw file
  count_bw <- function(bw_path,chr_list) {
    bw <- BigWigFileList(bw_path)
    coverage <- import(bw[[1]], as = 'RleList')
    s_count <- data.frame()
    for (i in c(1:length(chr_list))) {
      temp_chr <- e_merge_by_chr[which(e_merge_by_chr$chr == chr_list[i]),]
      temp_ir <- IRanges(temp_chr$start, temp_chr$end, names = temp_chr$chr)
      temp_coverage <- coverage[[chr_list[i]]]
      count <- as.numeric(sum(Views(temp_coverage, ranges(temp_ir))))
      temp_out <- cbind(temp_chr$e_merge_name,count)
      colnames(temp_out) <- c("e_merge_name","count")
      rownames(temp_out) <- NULL
      s_count <- rbind(s_count,temp_out)
    }
    return(s_count)
  }

  # set mem.maxVSize to Inf otherwise will get memory limit error
  mem.maxVSize(vsize = Inf)
  chr_list <- unique(e_merge_by_chr$chr)
  s1_1_count <- count_bw(s1_r1_bw,chr_list)
  s1_2_count <- count_bw(s1_r2_bw,chr_list)
  s2_1_count <- count_bw(s2_r1_bw,chr_list)
  s2_2_count <- count_bw(s2_r2_bw,chr_list)

  # get final count matrix
  counts_merge_1 <- merge(s1_1_count,s1_2_count,by="e_merge_name")
  counts_merge_2 <- merge(s2_1_count,s2_1_count,by="e_merge_name")
  counts_merge <- merge(counts_merge_1,counts_merge_2,by="e_merge_name")
  count_matrix <- data.frame(counts_merge[,-1], row.names=counts_merge$e_merge_name)
  colnames(count_matrix) <- c("S1_r1","S1_r2","S2_r1","S2_r2")
  count_matrix$S1_r1 <- as.integer(count_matrix$S1_r1)
  count_matrix$S1_r2 <- as.integer(count_matrix$S1_r2)
  count_matrix$S2_r1 <- as.integer(count_matrix$S2_r1)
  count_matrix$S2_r2 <- as.integer(count_matrix$S2_r2)

  # make sample data for DESeq2
  sample_data <- data.frame(row.names = c("S1_r1","S1_r2",
                                          "S2_r1","S2_r2"),
                            condition = c("S1","S1","S2","S2"))

  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = sample_data,
                                design = ~ condition)

  # run DEseq2
  dds <- DESeq(dds)
  res <- results(dds)

  # save raw counts and normalized counts

  raw_count <- counts(dds)
  normalized_count <- counts(dds,normalized=TRUE)
  colnames(normalized_count) <- c("S1_r1_norm","S1_r2_norm","S2_r1_norm","S2_r2_norm")

  # shrinking LFC and MA plot
  resLFC <- lfcShrink(dds, coef="condition_S2_vs_S1", type="apeglm")

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

  # make output list
  out <- list()
  out$enhancer_deseq_result <- final_out_df
  out$size_factor <- size_factor
  return(out)
}
