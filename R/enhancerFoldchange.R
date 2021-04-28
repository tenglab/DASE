#' Calculate the fold change of enhancers
#'
#' A function to calculate the fold change of enhancers base on two conditions.
#'
#' @details
#' This function will extract the enhancers which in the range of SEs.
#' Then enhancers are merged if there gaps are within 500bps.
#' Feature count and DESeq2 are used to calculate the fold change of those enhancers
#'
#' @param e_df pooled macs2 output file *_peaks.narrowPeak of all samples.
#' Header needs to be "chr, start,end,name,score,strand,signalValue,pValue,qValue,peak"
#'
#' @param se_df merged SE metadata from SEfilter
#' @param s1_r1_bam path of sample 1 replicate 1 bam file
#' @param s1_r2_bam path of sample 1 replicate 2 bam file
#' @param s2_r1_bam path of sample 2 replicate 1 bam file
#' @param s2_r2_bam path of sample 2 replicate 2 bam file
#' @param s1_pair if sample 1 is paired-end (default=FALSE)
#' @param s2_pair if sample 2 is paired-end (default=FALSE)
#'
#' @return
#' A list of 2 datasets: enhancer foldchang results dataset and sizefactor used in DESeq2
#'
#' @import data.table
#' @import DESeq2
#' @import Rsubread
#' @import apeglm
#'
#' @export
#' @examples
#' foldchange_list <- enhancerFoldchange(pool_enhancer_df,se_meta,s1_r1_bam,s1_r2_bam,s2_r1_bam,s2_r2_bam)
#'

enhancerFoldchange <- function(e_df,se_df,
                               s1_pair=FALSE,s2_pair=FALSE,
                               s1_r1_bam,s1_r2_bam,s2_r1_bam,s2_r2_bam) {

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


  # feature count
  # create SAF format file
  fc_saf <- e_merge_by_chr[,c(5,1,2,3)]
  fc_saf$Strand <- rep(".",nrow(fc_saf))
  colnames(fc_saf) <- c("GeneID","Chr","Start","End","Strand")

  # count each bam files
  # check if sample 1 is  paired-end
  if (s1_pair == F) {
    s1_r1_fc <- featureCounts(files=s1_r1_bam,annot.ext=fc_saf,
                              isGTFAnnotationFile = "SAF")
    s1_r2_fc <- featureCounts(files=s1_r2_bam,annot.ext=fc_saf,
                                  isGTFAnnotationFile = "SAF")
  } else {
    s1_r1_fc <- featureCounts(files=s1_r1_bam,annot.ext=fc_saf,
                              isGTFAnnotationFile = "SAF",
                              isPairedEnd=TRUE)
    s1_r2_fc <- featureCounts(files=s1_r2_bam,annot.ext=fc_saf,
                              isGTFAnnotationFile = "SAF",
                              isPairedEnd=TRUE)
  }

  # check if smaple 2 is paired-end
  if (s2_pair == F) {
    s2_r1_fc <- featureCounts(files=s2_r1_bam,annot.ext=fc_saf,
                              isGTFAnnotationFile = "SAF")
    s2_r2_fc <- featureCounts(files=s2_r2_bam,annot.ext=fc_saf,
                                  isGTFAnnotationFile = "SAF")
  } else {
    s2_r1_fc <- featureCounts(files=s2_r1_bam,annot.ext=fc_saf,
                              isGTFAnnotationFile = "SAF",
                              isPairedEnd=TRUE)
    s2_r2_fc <- featureCounts(files=s2_r2_bam,annot.ext=fc_saf,
                              isGTFAnnotationFile = "SAF",
                              isPairedEnd=TRUE)
  }

  # DESeq2 fold-change
  # make count matrix
  count_matrix <- as.data.frame(cbind(s1_r1_fc$counts,s1_r2_fc$counts,
                                      s2_r1_fc$counts,s2_r2_fc$counts))
  colnames(count_matrix) <- c("S1_r1","S1_r2","S2_r1","S2_r2")

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
  out$lfc_shrink <- resLFC
  out$enhancer_deseq_result <- final_out_df
  out$size_factor <- size_factor
  return(out)
}
