#!/usr/bin/bash -l

#SBATCH --time=48:00:00
#SBATCH --mem=10g                                                                                  
#SBATCH --cpus-per-task=1




# needed for Athef because it didn't automatically source ~/.bashrc so couldn't conda activate anything
source /home/hsiehph/shared/bin/initialize_conda.sh

module purge

conda activate snakemake


snakemake -s repeatmasker.snake --jobname "{rulename}.{jobid}" --profile profile  -j 100 -k --rerun-incomplete --restart-times 1 --use-conda \
          --groups MaskContig=MaskContiggroup \
          SpecialMaskContig=SpecialMaskContiggroup \
          TRFMaskContig=TRFMaskContiggroup \
          MergeMaskerRuns=MergeMaskerRunsgroup \
          CombineMaskedFasta=CombineMaskedFastagroup \
          FixCoordinates=FixCoordinatesgroup \
          --group-components MaskContiggroup=16  \
           SpecialMaskContiggroup=16 \
           TRFMaskContiggroup=16 \
           MergeMaskerRunsgroup=16 \
           CombineMaskedFastagroup=16 \
           FixCoordinatesgroup=16
