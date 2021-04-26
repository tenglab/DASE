# SEprofiler
R package used to categorize and rank super-enhancer (SE) subtle differences in two conditions.

SEprofiler R package:

Version 0.1.0

Dependency:

DESeq2, GenomicRanges, Rsubread, apeglm, data.table, ggplot2, splines

## Install
```R
# Public package
devtools::install_github("https://github.com/tenglab/SEprofiler.git")

# Private with token
devtools::install_github("https://github.com/tenglab/SEprofiler.git",auth_toke="abc")

# load package
library(SEprofiler)
```

## Files need to pre-process before using package                                                               
1. merge ROSE *_peaks_Gateway_SuperEnhancers.bed output file                                                 
2. merge macs2 *_peaks.narrowPeak output output file                                                         
3. path to sorted bam files of each sample and replicates                                                                     

## Pre-process
### 1. ROSE SE output bed file
```R
se_sample_1_r1 <- read.table("input/BC1_1_peaks_Gateway_SuperEnhancers.bed",sep='\t', header =F)
se_sample_1_r2 <- read.table("input/BC1_2_peaks_Gateway_SuperEnhancers.bed",sep='\t', header =F)
se_sample_2_r1 <- read.table("input/BC3_1_peaks_Gateway_SuperEnhancers.bed",sep='\t', header =F)
se_sample_2_r2 <- read.table("input/BC3_2_peaks_Gateway_SuperEnhancers.bed",sep='\t', header =F)

# merge sample and change the column name exact as follow
pool_se_df <- rbind(se_sample_1_r1,se_sample_2_r1,se_sample_1_r2,se_sample_2_r2)
colnames(pool_se_df) <- c("CHROM","START","STOP","REGION_ID","Signal","STRAND")
```

### 2. macs2 output peak file
```R
enhancer_s1_r1 <- read.table("input/BC1_1_peaks.narrowPeak",sep='\t')
enhancer_s1_r2 <- read.table("input/BC1_2_peaks.narrowPeak",sep='\t')
enhancer_s2_r1 <- read.table("input/BC3_1_peaks.narrowPeak",sep='\t')
enhancer_s2_r2 <- read.table("input/BC3_2_peaks.narrowPeak",sep='\t')

# merge sample and change the column name exact as follow
pool_enhancer_df <- rbind(enhancer_s1_r1,enhancer_s1_r2,enhancer_s2_r1,enhancer_s2_r2)
colnames(pool_enhancer_df) <- c("chr", "start","end","name","score","strand", "signalValue","pValue","qValue","peak")
```

### 3. path to bam file
```R
s1_r1_bam_path <- "input/BC1_1_rmdup_sort.bam"
s1_r2_bam_path <- "input/BC1_2_rmdup_sort.bam"
s2_r1_bam_path <- "input/BC3_1_rmdup_sort.bam"
s2_r2_bam_path <- "input/BC3_2_rmdup_sort.bam"
```

### 4. super-enhancer blacklist file (optional)
```R
blacklist_df <- read.table("input/ENCFF356LFX_blacklist.bed",sep = '\t')
```

# Usage
## Input file                                                                                                 
 1. pool_se_df                                                                                                     
 2. pool_enhancer_df                                                                                          
 3. s1_r1_bam_path, s1_r2_bam_path, s2_r1_bam_path, s2_r2_bam_path                                             

## Examples
```R
# default: no blacklist file, no permutation, default cutoff=c(-1.5,1.5), both samples are single end
se_test_out <- SEprofile(se_in = pool_se_df, e_df = pool_enhancer_df, 
                         s1_r1_bam = s1_r1_bam_path, s1_r2_bam = s1_r2_bam_path, 
                         s2_r1_bam = s2_r1_bam_path, s2_r2_bam = s2_r2_bam_path)

# with blacklist, no permutation, personal defined cutoff, both samples are paired-end
se_test_out <- SEprofile(se_in = pool_se_df, e_df = pool_enhancer_df, 
                         bl_file = blacklist_df, has_bl_file = TRUE,
                         cutoff_v = c(-1.7,1.7),
                         s1_pair = T, s2_pair = T,
                         s1_r1_bam = s1_r1_bam_path, s1_r2_bam = s1_r2_bam_path,
                         s2_r1_bam = s2_r1_bam_path, s2_r2_bam = s2_r2_bam_path)
                           
# with blacklist, with permutation (cutoff will defined by permutation), sample 1 is paired-end
se_test_out <- SEprofile(se_in = pool_se_df, e_df = pool_enhancer_df, 
                         bl_file = blacklist_df, has_bl_file = T, permut = T,
                         s1_pair = T, s2_pair = F,
                         s1_r1_bam = s1_r1_bam_path, s1_r2_bam = s1_r2_bam_path,
                         s2_r1_bam = s2_r1_bam_path, s2_r2_bam = s2_r2_bam_path)
                         
# save pattern plots to pdf
dir.create("output/patterns")
pattern_list <- se_profile_out$pattern_plot

for (i in c(1:length(pattern_list))) {
  # create out plot name
  pdf_name <- paste(gsub(" pattern", "",pattern_list[[i]]$labels$title), ".pdf", sep="")
  
  pdf_full_path <- paste("output/patterns",pdf_name,sep="/")
  
  # save to pdf
  pdf(file = pdf_full_path,width = 8, height = 4)
    print(pattern_list[[i]])
  dev.off()
}

```
## Output list
 1. se_test_out$cate_rank: final SE categories and ranking
 2. se_test_out$box_plot: boxplot of super-enhancer categories
 3. se_test_out$permut_plot: permutation density plot if permut=T, return NA if permut=F
 4. se_test_out$fc_cutoff: fold change cutoffs
 5. se_test_out$e_fit_fc: counts fold change and fitted fold change of enhancer within SEs 
 6. se_test_out$pattern_plot: list contain fitted line pattern plot of each SE
 7. se_test_out$se_segments: SE segment percentage
