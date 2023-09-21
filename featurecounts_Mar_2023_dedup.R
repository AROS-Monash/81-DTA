# author: Shani Amarasinghe
# affiliation: Rosello-Diez lab, ARMI, Monash
# date: 7 July 2023
# analysing Trev-seq bulk RNA-Seq data for Chee - M3
# generating count matrix for the aligned reads using Rsubsread::featureCounts()
# for the new 2023 alignment - and dedupl using umitools

# load R and libraries
library(Rsubread)
library(ggplot2)

# load I/O locations
bam_dir   <-"/fs03/ag36/Shani/Trev-Seq/STAR_aligned_Mar_2023_dedup/"
bam.files <- list.files(bam_dir, pattern = ".bam$", full.names =TRUE)
GTF       <- "/fs04/ag36/Shani/references/Mus_musculus.GRCm39.109.chr.gtf"

fc <- featureCounts(bam.files,

    # annotation
    annot.ext = GTF,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    
    # level of summarization
    useMetaFeatures = TRUE,
    
    # overlap between reads and features
    allowMultiOverlap = TRUE,
    minOverlap = 1,
    fracOverlap = 0,
    fracOverlapFeature = 0,
    largestOverlap = FALSE,
    nonOverlap = NULL,
    nonOverlapFeature = NULL,

    # Read shift, extension and reduction
    readShiftType = "upstream",
    readShiftSize = 0,
    readExtension5 = 0,
    readExtension3 = 0,
    read2pos = NULL,
    
    # multi-mapping reads - this was tested for both TRUE and FALSE
    countMultiMappingReads = TRUE,

    # fractional counting
    fraction = FALSE,

    # long reads
    isLongRead = FALSE,

    # read filtering
    minMQS = 0,
    splitOnly = FALSE,
    nonSplitOnly = FALSE,
    primaryOnly = FALSE,
    ignoreDup = FALSE,
    
    # strandness
    strandSpecific = 0,
    
    # exon-exon junctions
    juncCounts = FALSE,
    genome = NULL,
    
    # parameters specific to paired end reads
    isPairedEnd = FALSE,
    countReadPairs = TRUE,
    requireBothEndsMapped = FALSE,
    checkFragLength = FALSE,
    minFragLength = 50,
    maxFragLength = 600,
    countChimericFragments = TRUE,    
    autosort = TRUE,
    
    # number of CPU threads
    nthreads = 12,

    # read group
    byReadGroup = FALSE,

    # report assignment result for each read
    reportReads = NULL,
    reportReadsPath = NULL,

    # miscellaneous
    maxMOp = 10,
    tmpDir = ".",
    verbose = TRUE)

# saving multimaping allowed
write.table(fc$stat, "/fs03/ag36/Shani/Trev-Seq/Counts_matrix/trev-seq-mm-countmatrix-stats_Mar2023_dedup.tab", quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)
write.table(fc$counts, "/fs03/ag36/Shani/Trev-Seq/Counts_matrix/trev-seq-mm-countmatrix_Mar2023_dedup.tab", quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)
saveRDS(fc, "/fs03/ag36/Shani/Trev-Seq/Counts_matrix/trev-seq-mm-countmatrix_Mar2023_dedup.RDS")

fc.nomm <- featureCounts(bam.files,

    # annotation
    annot.ext = GTF,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    
    # level of summarization
    useMetaFeatures = TRUE,
    
    # overlap between reads and features
    allowMultiOverlap = TRUE,
    minOverlap = 1,
    fracOverlap = 0,
    fracOverlapFeature = 0,
    largestOverlap = FALSE,
    nonOverlap = NULL,
    nonOverlapFeature = NULL,

    # Read shift, extension and reduction
    readShiftType = "upstream",
    readShiftSize = 0,
    readExtension5 = 0,
    readExtension3 = 0,
    read2pos = NULL,
    
    # multi-mapping reads - this was tested for both TRUE and FALSE
    countMultiMappingReads = FALSE,

    # fractional counting
    fraction = FALSE,

    # long reads
    isLongRead = FALSE,

    # read filtering
    minMQS = 0,
    splitOnly = FALSE,
    nonSplitOnly = FALSE,
    primaryOnly = FALSE,
    ignoreDup = FALSE,
    
    # strandness
    strandSpecific = 0,
    
    # exon-exon junctions
    juncCounts = FALSE,
    genome = NULL,
    
    # parameters specific to paired end reads
    isPairedEnd = FALSE,
    countReadPairs = TRUE,
    requireBothEndsMapped = FALSE,
    checkFragLength = FALSE,
    minFragLength = 50,
    maxFragLength = 600,
    countChimericFragments = TRUE,    
    autosort = TRUE,
    
    # number of CPU threads
    nthreads = 12,

    # read group
    byReadGroup = FALSE,

    # report assignment result for each read
    reportReads = NULL,
    reportReadsPath = NULL,

    # miscellaneous
    maxMOp = 10,
    tmpDir = ".",
    verbose = TRUE)


# saving multimapping not allowed
write.table(fc.nomm$stat, "/fs03/ag36/Shani/Trev-Seq/Counts_matrix/trev-seq-nomm-countmatrix-stats_Mar2023_dedup.tab")
write.table(fc.nomm$counts, "/fs03/ag36/Shani/Trev-Seq/Counts_matrix/trev-seq-nomm-countmatrix_Mar2023_dedup.tab")
saveRDS(fc.nomm, "/fs03/ag36/Shani/Trev-Seq/Counts_matrix/trev-seq-nomm-countmatrix_Mar2023_dedup.RDS")

