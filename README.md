# DASE
R package used to categorize and rank super-enhancer (SE) subtle differences in two conditions.

DASE R package:

Version 0.1.0

Dependency:

DESeq2, GenomicRanges, Rsubread, apeglm, data.table, ggplot2, splines

## Install
```R
# Public package
devtools::install_github("https://github.com/tenglab/DASE.git")

# Private with token
devtools::install_github("https://github.com/tenglab/DASE.git",auth_toke="your token")

# load package
library(DASE)
```

## Files need to pre-process before using package                                                               
1. merge ROSE *_peaks_Gateway_SuperEnhancers.bed output file                                                 
2. merge macs2 *_peaks.narrowPeak output output file                                                         
3. path to sorted bam files or bigwig files of each sample and replicates                                                                     


# Usage
## Input file                                                                                                 
 1. pool_se_df                                                                                                     
 2. pool_enhancer_df                                                                                          
 3. s1_r1_bam_path, s1_r2_bam_path, s2_r1_bam_path, s2_r2_bam_path                                             

## Examples
```R
# default: no blacklist file, no permutation, default cutoff=c(-1.5,1.5), both samples are single end
se_test_out <- DASE(se_in = pool_se_df, e_df = pool_enhancer_df, 
                         s1_r1_bam = s1_r1_bam_path, s1_r2_bam = s1_r2_bam_path, 
                         s2_r1_bam = s2_r1_bam_path, s2_r2_bam = s2_r2_bam_path)

# with blacklist and customize blacklist range, no permutation, personal defined cutoff, both samples are paired-end
se_test_out <- DASE(se_in = pool_se_df, e_df = pool_enhancer_df, 
                         bl_file = blacklist_df, has_bl_file = TRUE, custom_range = c("chr1:123-12345","chr12:12345-1234567"),
                         cutoff_v = c(-1.7,1.7),
                         s1_pair = T, s2_pair = T,
                         s1_r1_bam = s1_r1_bam_path, s1_r2_bam = s1_r2_bam_path,
                         s2_r1_bam = s2_r1_bam_path, s2_r2_bam = s2_r2_bam_path)
                           
# with blacklist, with permutation (cutoff will defined by permutation),and bw coverage
se_test_out <- DASE(se_in = pool_se_df, e_df = pool_enhancer_df, 
                         bl_file = blacklist_df, has_bl_file = T, permut = T, data_type="bw"
                         s1_r1_bam = s1_r1_bw_path, s1_r2_bam = s1_r2_bw_path,
                         s2_r1_bam = s2_r1_bw_path, s2_r2_bam = s2_r2_bw_path)
                         
# save pattern plots to pdf
dir.create("output/patterns")
pattern_list <- se_test_out$pattern_list

for (i in c(1:length(pattern_list))) {
  # create out plot name
  cat_dir <- se_profile_out$se_category$category[which(se_test_out$se_category$se_merge_name==
                                                         pattern_list[[i]]$labels$title)]
  pdf_name <- paste(gsub(" pattern", "",pattern_list[[i]]$labels$title),".pdf",sep="")
  pdf_full_path <- paste0("output_new/patterns/",cat_dir,"/",pdf_name)

  # save to pdf
  pdf(file = pdf_full_path,width = 8, height = 4)
  print(pattern_list[[i]])
  dev.off()
}

```
## Output list
 1. se_test_out$se_category: final SE categories and ranking
 2. se_test_out$boxplot: boxplot of super-enhancer categories
 3. se_test_out$density_plot: permutation density plot if permut=T, return NA if permut=F
 4. se_test_out$cutoff: fold change cutoffs
 5. se_test_out$se_fit: counts fold change and fitted fold change of enhancer within SEs 
 6. se_test_out$pattern_list: list contain fitted line pattern plot of each SE
 7. se_test_out$se_seg: SE segment percentage
 8. se_test_out$se_meta: merged SE bed
 9. se_test_out$se_deseq_out: constitute enhancer DEseq2 out table
