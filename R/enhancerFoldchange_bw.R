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
#'
#' @param se_df merged SE metadata from SEfilter
#' @param condition_1 path of bam files for condition 1 replicates, separated by ",".
#'                    Example: condition_1="path_to_bw/c1_r1.bw,path_to_bw/c1_r2.bw,...".
#' @param condition_2 path of bam files for condition 1 replicates, separated by ",".
#'                    Example: condition_2="path_to_bw/c2_r1.bw,path_to_bw/c2_r2.bw,...".
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
#' foldchange_list <- enhancerFoldchange(pool_enhancer_df,se_meta,
#'                        condition_1="path_to_bw/c1_r1.bw,path_to_bam/c1_r2.bw,...",
#'                        condition_2="path_to_bw/c2_r1.bw,path_to_bam/c2_r2.bw,...")
#'

enhancerFoldchange_bw <- function(e_df,se_df,condition_1,condition_2) {

  # merge enhancers with gaps less than 500bps
  ir <- IRanges(e_df$V2,
                e_df$V3,
                names = e_df$V1)

  e_merge_by_chr <- rbindlist(lapply(split(ir, names(ir)),
                                     function(x) as.data.table(reduce(x,min.gapwidth = 500))),
                              idcol = "chr")

  colnames(e_merge_by_chr) <- c("chr","start","end","width")

  # create merged enhancer names
  e_merge_by_chr$e_merge_name <- paste0(e_merge_by_chr$chr,"_",e_merge_by_chr$start,"_",e_merge_by_chr$end)

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


  # get each replicates and if they are paired-end reads
  c1_bw <- unlist(strsplit(gsub(" ","",condition_1),","))
  c2_bw <- unlist(strsplit(gsub(" ","",condition_2),","))

  # count each replicates
  chr_list <- unique(e_merge_by_chr$chr)

  # check number of replicates >= 2
  if (length(c1_bw)>=2 & length(c2_bw)>=2) {
    # initial count matrix
    count_matrix_c1 <- data.frame(matrix(nrow = nrow(e_merge_by_chr),ncol = length(c1_bw)))
    count_matrix_c2 <- data.frame(matrix(nrow = nrow(e_merge_by_chr),ncol = length(c2_bw)))
    row.names(count_matrix_c1) <- e_merge_by_chr$e_merge_name
    row.names(count_matrix_c2) <- e_merge_by_chr$e_merge_name

    # count condition 1
    for (c1_i in 1:length(c1_bw)) {
      c1_fc_tmp <- paste0("C1_",c1_i)
      assign(c1_fc_tmp,count_bw(c1_bw[c1_i],chr_list))
      count_matrix_c1[,c1_i] <- as.numeric(get(c1_fc_tmp)[,2])
      colnames(count_matrix_c1)[c1_i] <- c1_fc_tmp
    }

    # count condition 2
    for (c2_i in 1:length(c2_bw)) {
      c2_fc_tmp <- paste0("C2_",c2_i)
      assign(c2_fc_tmp,count_bw(c2_bw[c2_i],chr_list))
      count_matrix_c2[,c2_i] <- as.numeric(get(c2_fc_tmp)[,2])
      colnames(count_matrix_c2)[c2_i] <- c2_fc_tmp
    }

  } else {
    print("DASE need more than 2 replicates in each condition.")
  }

  # DESeq2 fold-change for enhancer in SE only
  # make count matrix
  count_matrix <- cbind(count_matrix_c1,count_matrix_c2)

  # get final count matrix
  # remove enhancers with 0 value for all sample
  not_zero_in_se <- apply(count_matrix, 1, function(row) sum(row) != 0)
  count_matrix <- count_matrix[not_zero_in_se,]

  # make sample data for DESeq2
  sample_data <- data.frame(row.names = colnames(count_matrix),
                            condition = c(rep("C1",length(c1_bw)),rep("C2",length(c2_bw))))

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

