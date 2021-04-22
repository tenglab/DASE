# SEprofiler
R package used to category and ranking super-enhancer (SE) subtle differences in two conditions.

SEprofiler R package:

Version 0.1.0

github: https://github.com/tenglab/SEprofiler

Dependency:

DESeq2, GenomicRanges, Rsubread, apeglm, data.table, ggplot2, splines


# Files need to pre-process before using package                                                               
1. merge ROSE *_peaks_Gateway_SuperEnhancers.bed output file                                                 
2. merge macs2 *_peaks.narrowPeak output output file                                                         
3. path to sorted bam files of each sample and replicates                                                    
                                                                                                              
# options:                                                                                                     
super-enhancer blacklist file                                                                                

# pre-process example
# rose SE output bed file
main_sample_1_r1 &lt;- read.table("input/BC1_1_peaks_Gateway_SuperEnhancers.bed",sep='\t', header =F)
main_sample_1_r2 <- read.table("input/BC1_2_peaks_Gateway_SuperEnhancers.bed",sep='\t', header =F)
main_sample_2_r1 <- read.table("input/BC3_1_peaks_Gateway_SuperEnhancers.bed",sep='\t', header =F)
main_sample_2_r2 <- read.table("input/BC3_2_peaks_Gateway_SuperEnhancers.bed",sep='\t', header =F)

# merge sample and change the column name exact as follow
se_df <- rbind(main_sample_1_r1,main_sample_2_r1,main_sample_1_r2,main_sample_2_r2)
colnames(se_df) <- c("CHROM","START","STOP","REGION_ID","Signal","STRAND")

# super-enhancer blacklist file (optional)
blacklist_df <- read.table("../../ENCFF356LFX_blacklist.bed",sep = '\t')

# macs2 output peak file
enhancer_s1_r1 <- read.table("input/BC1_1_peaks.narrowPeak",sep='\t')
enhancer_s1_r2 <- read.table("input/BC1_2_peaks.narrowPeak",sep='\t')
enhancer_s2_r1 <- read.table("input/BC3_1_peaks.narrowPeak",sep='\t')
enhancer_s2_r2 <- read.table("input/BC3_2_peaks.narrowPeak",sep='\t')

# merge sample and change the column name exact as follow
pool_enhancer_df <- rbind(enhancer_s1_r1,enhancer_s1_r2,enhancer_s2_r1,enhancer_s2_r2)
colnames(pool_enhancer_df) <- c("chr", "start","end","name","score","strand", "signalValue","pValue","qValue","peak")


# path to bam file
s1_r1_bam_path <- "~/Projects/super_enhancer/SE_paper_data/PEL_GSE136090/input/BC1_1_rmdup_sort.bam"
s1_r2_bam_path  <- "~/Projects/super_enhancer/SE_paper_data/PEL_GSE136090/input/BC1_2_rmdup_sort.bam"
s2_r1_bam_path <- "~/Projects/super_enhancer/SE_paper_data/PEL_GSE136090/input/BC3_1_rmdup_sort.bam"
s2_r2_bam_path <- "~/Projects/super_enhancer/SE_paper_data/PEL_GSE136090/input/BC3_2_rmdup_sort.bam"


# Input file:                                                                                                  
 1. se_df                                                                                                     
 2. pool_enhancer_df                                                                                          
 3. s1_r1_bam_path, s1_r2_bam_path, s2_r1_bam_path, s2_r2_bam_path                                             
                                                                                                              
# options:                                                                                                     
 1. blacklist_df                                                                                              
 2. if using permutation (Default fold change cutoff is (-1.5,1.5) if not using permutation)


# Using example:
# default: no blacklist file, permutation 10 times

se_test_out <- SEmain(se_in = se_df, e_df = pool_enhancer_df, 
                      s1_r1_bam = s1_r1_bam_path, s1_r2_bam = s1_r2_bam_path,
                      s2_r1_bam = s2_r1_bam_path, s2_r2_bam = s2_r2_bam_path)

# with blacklist file, no permutation
se_test_out <- SEmain(se_in = se_df, e_df = pool_enhancer_df, 
                      bl_file = blacklist_df,has_bl_file = TRUE, permut = FALSE,
                      s1_r1_bam = s1_r1_bam_path, s1_r2_bam = s1_r2_bam_path,
                      s2_r1_bam = s2_r1_bam_path, s2_r2_bam = s2_r2_bam_path)
