#!/bin/bash
#SBATCH --job-name=limma-dedup    # Job name
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=shani.amarasinghe@monash.edu     # Where to send mail	
#SBATCH --ntasks=18                    # Run on a single CPU
#SBATCH --mem=160gb                     # Job memory request
#SBATCH --time=2-20:23:00               # Time limit hrs:min:sec
#SBATCH --output=/fs04/ag36/Shani/Trev-Seq/limma-2023mar-dedup-slurm-%A_%a.out # Standard output and error log

pwd; hostname; date

module load R/4.2.2-mkl

Rscript --vanilla /fs04/ag36/Shani/Trev-Seq/scripts/edgeR-limma-genelevel-stage-CC-Mar-2023-dedup.R