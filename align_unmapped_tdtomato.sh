#!/bin/bash
#SBATCH --job-name=tdtomato_star
#SBATCH --account=ag36
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=11
#SBATCH --mem-per-cpu=5555
#SBATCH --cpus-per-task=1
#SBATCH --output=/fs04/ag36/Shani/Trev-Seq/scripts/tdtomato-star-align-%A_%a.out

# author: Shani Amarasinghe
# affiliation: Rosello-Diez lab, ARMI, Monash
# date: 16 Mar 2023
# Align the unmapped reads from trev-seq dataset to tdtomato and see whether they map
# align using STAR

# Define I/O -----------------------
BAM_DIR="/fs03/ag36/Shani/Trev-Seq/STAR_aligned_2023"
DATA_DIR="/fs03/ag36/Shani/Trev-Seq/Star_unmapped"
OUTPUT_DIR="/fs03/ag36/Shani/Trev-Seq/Star_unmapped/star_aligned"
TMP_DIR="/fs03/ag36/Shani/Trev-Seq/Star_unmapped/star_aligned/tmp"
GENOME="/fs03/ag36/Shani/Trev-Seq/Star_unmapped/tdtomato_ref"

mkdir $OUTPUT_DIR
mkdir $GENOME

# STAR -------------------------

# load modules for STAR again (gcc is different to what is needed by umi-tools) 
module unload gcc
module load star/2.7.9a
module load samtools/1.9

# extract the unmapped reads from the BAM files and converting them to fastq ----------

# find "${BAM_DIR}" -name '*_Aligned.sortedByCoord.out.bam' -print0 | while IFS= read -r -d '' BAM
# do
#     prefix=${BAM##*/}; prefix=${prefix%%Aligned.sortedByCoord.out.bam};
#     echo -ne " \n ################# extracting from $BAM ######### \n\n";
#     samtools view -b -f 4 "${BAM}" > ${DATA_DIR}/unmapped_"${prefix}".bam
#     samtools fastq -F 256 ${DATA_DIR}/unmapped_"${prefix}".bam > ${DATA_DIR}/unmapped_"${prefix}".fastq
# done

# creating STAR genome index -----------
STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir ${GENOME} \
--genomeSAindexNbases 4 \
--genomeFastaFiles /fs03/ag36/Shani/Trev-Seq/Star_unmapped/tdt.fasta

# # ------ started mapping with more customisation to adjust for smaller genomes and no SJs
# #--sjdbOverhang 99 
# #--sjdbGTFfile /fs03/ag36/Shani/Ehsan_knee_scRNASeq/STAR_aligned_GOI/transgenes.gtf \


# star mapping 

find "${DATA_DIR}" -name '*.fastq' -print0 | while IFS= read -r -d '' R1
do
    # extract the first part of the read name for output name
    prefix=${R1##*/}; prefix=${prefix%%.fastq*};
    echo -ne " \n ################# aligning $R1 ######### \n\n";
    # start mapping
    STAR \
    --genomeDir ${GENOME} \
    --outWigType bedGraph \
    --readFilesIn ${R1} \
    --genomeChrBinNbits 7.55 \
    --runThreadN 18 \
    --alignMatesGapMax 10 \
    --alignIntronMax 1 \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 \
    --outFileNamePrefix ${OUTPUT_DIR}/${prefix}_remapped_ \
    --outFilterMultimapNmax 0 \
    --outFilterMismatchNmax 2 \
    --outSAMunmapped Within \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattrRGline ID:${prefix} DS:  PE:Illumina PU:10XsnRNA-Seq SM:${prefix} \
    --outTmpDir ${TMP_DIR};
    echo -ne " \n################# complete aligning $R1 ######### \n\n";
    rm -r ${TMP_DIR};
done;