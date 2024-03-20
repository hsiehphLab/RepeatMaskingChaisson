#!/usr/bin/bash -l

#SBATCH --time=48:00:00
#SBATCH --mem=10g                                                                                  
#SBATCH --cpus-per-task=1




# needed for Athef because it didn't automatically source ~/.bashrc so couldn't conda activate anything
source /home/hsiehph/shared/bin/initialize_conda.sh

module purge

conda activate snakemake


snakemake -s RepeatMaskGenome.epp.snakefile --jobname "{rulename}.{jobid}" --profile profile  -j 500 -k --rerun-incomplete --restart-times 1 --use-conda
