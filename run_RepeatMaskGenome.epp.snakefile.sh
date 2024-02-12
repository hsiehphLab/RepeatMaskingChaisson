#!/usr/bin/bash -l

# needed for Athef because it didn't automatically source ~/.bashrc so couldn't conda activate anything
source ~/.bashrc

module purge

conda activate snakemake


snakemake -s RepeatMaskGenome.epp.snakefile --jobname "{rulename}.{jobid}" --profile profile  -j 500 -k --rerun-incomplete --restart-times 1 --use-conda