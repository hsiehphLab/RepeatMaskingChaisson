#!/usr/bin/env python
import sys
asm=sys.argv[1]
hap=sys.argv[2]
outFileName=sys.argv[3]
fileStr="""{{
   "asm" : "{}",
   "hap": "{}",
   "repeat_library": "human",
   "grid_large": "sbatch -c 16 --mem=48G --time=24:00:00 --partition=qcb --account=mchaisso_100", 
   "grid_medium": "sbatch -c 8 --mem=16G --time=24:00:00 --partition=qcb --account=mchaisso_100", 
   "grid_small": "sbatch -c 1 --mem=8G --time=4:00:00 --partition=qcb  --account=mchaisso_100", 
   "grid_blat": "sbatch -c 1 --mem=8G --time=24:00:00 --partition=qcb --account=mchaisso_100", 
   "grid_sedef": "sbatch -c 16 --mem=48G --time=24:00:00 --partition=qcb --account=mchaisso_100  --constraint=xeon-6130", 
   "grid_repeatmasker" : "sbatch -c 4 --mem=8G --time=24:00:00 --partition=qcb --account=mchaisso_100",
   "t2t_repeat_library" : "/project/mchaisso_100/projects/HPRC/Assemblies_v2/RepeatMaskerDatabase/final_consensi_gap_nohsat_teucer.embl.txt.fasta"
}}
""".format(asm,hap)

outFile=open(outFileName, 'w')
outFile.write(fileStr)
outFile.close()
            
