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

# Usage Examples
## Input file                                                                                                 
 1. SE peaks of all samples and replicates: test_data/se_peaks.bed
 
Header must be the same as follow
| CHROM | START     | STOP      | REGION_ID                  | Signal | STRAND |
| ----- |-----------| --------  | -------------------------- | ------ | ------ |
| chr14 | 34355128  | 34380484  | 5_Peak_2778_lociStitched   | 496    | .      |
| chr6  | 151869032 | 151917419 | 13_Peak_66874_lociStitched | 675    | .      |
| chr3  | 149297239 | 149317229 | 4_Peak_182_lociStitched    | 340    | .      |


 2. enhancer peaks of all samples and replicates: test_data/enhancer_peaks.bed  
 3. path of bam files: s1_r1_bam_path, s1_r2_bam_path, s2_r1_bam_path, s2_r2_bam_path
 4. SE blacklist (optional): test_data/ENCFF356LFX_blacklist.bed                                                                                  
                                             

```R
# default: no blacklist file, permutation 10 times, both samples are single end
se_test_out <- DASE(se_in = se_peaks.bed, e_df = enhancer_peaks.bed, 
                    s1_r1_bam = s1_r1_bam_path, s1_r2_bam = s1_r2_bam_path, 
                    s2_r1_bam = s2_r1_bam_path, s2_r2_bam = s2_r2_bam_path)

# with blacklist and customize blacklist range
se_test_out <- DASE(se_in = se_peaks.bed, e_df = enhancer_peaks.bed, 
                         bl_file = blacklist_df, has_bl_file = TRUE, custom_range = c("chr1:123-12345","chr12:12345-1234567"),
                         s1_r1_bam = s1_r1_bam_path, s1_r2_bam = s1_r2_bam_path,
                         s2_r1_bam = s2_r1_bam_path, s2_r2_bam = s2_r2_bam_path)
                           
# coverage is bw instead of bam
se_test_out <- DASE(se_in = se_peaks.bed, e_df = enhancer_peaks.bed, 
                         bl_file = blacklist_df, has_bl_file = T, permut = T, data_type="bw"
                         s1_r1_bam = s1_r1_bw_path, s1_r2_bam = s1_r2_bw_path,
                         s2_r1_bam = s2_r1_bw_path, s2_r2_bam = s2_r2_bw_path)
                         

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
