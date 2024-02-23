import os
import tempfile
import subprocess
import os.path


# Config
configfile: "repeatmask.json"

assembly=config["input"]
asmFai=assembly + ".fai"

if not os.path.exists( asmFai ):
	szCommand = "module load samtools/1.9 && samtools faidx " + assembly
	print( "about to execute: " + szCommand )
	subprocess.call( szCommand, shell = True )


asmFaiFile=open(asmFai)
faiLines=asmFaiFile.readlines()

szMaskedAssembly=config["output"]
szCompressedMaskedAssembly = szMaskedAssembly + ".gz"
szLocationsOfRepeats=config["output"] + ".out"

contigs=[l.split()[0] for l in faiLines]
lengths=[int(l.split()[1]) for l in faiLines]
strIdx=[str(i) for i in range(0,len(contigs))]


# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

SplitSize=2400000
SplitOverlap=20000

splitRegions=[]
splitByContig={}
contigNames=[]
#for i in range(0,len(contigs)):
for i in range(0,len(contigs)):

    name=contigs[i]
    contigNames.append(name)
    seqLen = lengths[i]
    splitByContig[name] = []
    nSeq = int(seqLen/SplitSize)
        
    if seqLen % SplitSize > 0:
        nSeq+=1
        
    for idx in range(0,nSeq):
        if idx == 0:
            start=0
            ovp=0
        else:
            start=idx*SplitSize-SplitOverlap
            ovp=SplitOverlap
        end=min((idx+1)*SplitSize, seqLen)

#        seqRec.id=seqRec.id.replace("/", "_").replace("|","_")
        coords=str(start) + "_" + str(end) + "_" + str(ovp)+ "_" + str(idx)
        rgn=name + "_" +coords

        splitByContig[name].append(coords)
        splitRegions.append( rgn)
sys.stderr.write("Processing {} regions\n".format(str(len(splitRegions))))






rule all:
	input:
		split=expand("split/to_mask.{region}.fasta",region=splitRegions),
		masked=expand("masked/to_mask.{region}.fasta.masked",region=splitRegions),
		t2t=expand("t2t/to_mask.{region}.fasta.masked",region=splitRegions),
		comb=expand("comb/to_mask.{region}.fasta.masked",region=splitRegions),
		mask=szCompressedMaskedAssembly,
		maskedGenomeOut=szLocationsOfRepeats
	shell:"""
		./cleanUp.py
"""


rule IndexGenome:
    input:
        asm=assembly
    output:
        fai=asmFai
    resources:
        mem=8,
        threads=1,
        hrs=4
    shell:"""
samtools faidx {input.asm}
"""

rule SplitGenome:
    input:
        asm=assembly,
        fai=asmFai
    output:
        split=expand("split/to_mask.{region}.fasta", region=splitRegions)
#    "grid_medium": "sbatch -c 8 --mem=16G --time=24:00:00 --partition=qcb --account=mchaisso_100", 
    resources:
        mem=16,
        threads=1,
        hrs=24
    params:
		#        grid_opts=config["grid_medium"],
        sd=SD
    resources:
        load=1
    shell:"""
mkdir -p split
module load python3/3.9.3_anaconda2021.11_mamba && {params.sd}/DivideFasta.py {input.asm} {SplitSize} {SplitOverlap} split/to_mask
"""

rule MaskContig:
    input:
        contig="split/to_mask.{index}.fasta"
    output:
        mask="masked/to_mask.{index}.fasta.masked",
        maskOut="masked/to_mask.{index}.fasta.out"
    conda:
        "rmsk"
    resources:
        mem=8,
        threads=1,
        hrs=24
    params:
        repeatLibrary=config["repeat_library"],
        sd=SD        
    shell:"""
mkdir -p masked
TEMP="$TMPDIR/$$_$RANDOM/"
mkdir -p $TEMP
cp \"{input.contig}\" \"$TEMP/to_mask.{wildcards.index}.fasta\" && \
pushd $TEMP && 
RepeatMasker -species "Homo sapiens" -pa 8  -s -xsmall \"to_mask.{wildcards.index}.fasta\" && \
popd && \
if [ ! -e $TEMP/to_mask.\"{wildcards.index}\".fasta.masked ]; then
  cp split/to_mask.\"{wildcards.index}\".fasta masked/to_mask.\"{wildcards.index}\".fasta.masked
  cp {params.sd}/repeat_masker.out.header masked/to_mask.\"{wildcards.index}\".fasta.out
else
  cp $TEMP/to_mask.\"{wildcards.index}\".fasta.* masked/ || true
fi
#rm -rf $TEMP
"""

rule SpecialMaskContig:
    input:
        mask="masked/to_mask.{index}.fasta.masked",
    output:
        tt="t2t/to_mask.{index}.fasta.masked",
        ttOut="t2t/to_mask.{index}.fasta.out",
    conda:
        "rmsk"
    resources:
        mem = 8,		
        threads = 1,
        hrs = 4
    params:
        repeatLibrary=config["t2t_repeat_library"],
        sd=SD
    shell:"""
if [ {params.repeatLibrary} != "na" ]
then
  mkdir -p t2t
  TEMP="$TMPDIR/$$_$RANDOM/"
  mkdir -p $TEMP
  # 
  # Copy input file to temp dir that should have fast IO
  #
  export LD_LIBRARY_PATH=/panfs/jay/groups/7/hsiehph/gordo893/packages/htslib && {params.sd}/hardmask {input.mask} $TEMP/to_mask.{wildcards.index}.fasta && \
  pushd $TEMP &&  \
  RepeatMasker -nolow -lib {params.sd}/T2TLib/final_consensi_gap_nohsat_teucer.embl.txt.fasta  -pa 8 -s -xsmall \"to_mask.{wildcards.index}.fasta\" && \
  popd && \

  if [ ! -e $TEMP/to_mask.\"{wildcards.index}\".fasta.masked ]; then
    ls -l masked/to_mask.\"{wildcards.index}\".fasta.masked
    cp masked/to_mask.\"{wildcards.index}\".fasta.masked t2t/to_mask.\"{wildcards.index}\".fasta.masked
    head -3 masked/to_mask.\"{wildcards.index}\".fasta.out > t2t/to_mask.\"{wildcards.index}\".fasta.out
  else
    cp $TEMP/to_mask.\"{wildcards.index}\".fasta.* t2t/ || true
  fi
  rm -rf $TEMP
else
  echo "*** rule SpecialMaskContig made no change (skipped) ***" > t2t/to_mask.\"{wildcards.index}\".fasta.out
  cp {input.mask} t2t/to_mask.\"{wildcards.index}\".fasta.masked >> t2t/to_mask.\"{wildcards.index}\".fasta.out || true
fi
"""
    
rule TRFMaskContig:
    input:
        orig="split/to_mask.{index}.fasta",
    output:
        trf="trf/to_mask.{index}.fasta.trf"
    conda:
        "base"
    resources:
        mem = 8,		
        threads = 1,
        hrs = 4
    params:
        sd=SD
    shell:"""
mkdir -p trf    
/panfs/jay/groups/7/hsiehph/shared/software/packages/trf-mod/from_Mark/TRF-mod/trf-mod -p 500 -l 20000 {input.orig} > {output.trf}.bed
# fix for /panfs/jay/groups/7/hsiehph/gordo893/pipelines/mark_chaisson_repeatmasker/RepeatMasking/bemask:     error while loading shared libraries: libhts.so.3: cannot open shared object file: No such file or directory
# fix for /usr/bin/bash: line 3: LD_LIBRARY_PATH: unbound variable
# and ./bemask: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found (required by ./bemask)
module load gcc/7.2.0 && export LD_LIBRARY_PATH=/panfs/jay/groups/7/hsiehph/gordo893/packages/htslib:$LD_LIBRARY_PATH && {params.sd}/bemask {input.orig} {output.trf}.bed {output.trf}
"""

rule MergeMaskerRuns:
    input:
        humLib="masked/to_mask.{index}.fasta.masked",
        humLibOut="masked/to_mask.{index}.fasta.out",
        t2tLib="t2t/to_mask.{index}.fasta.masked",
        t2tLibOut="t2t/to_mask.{index}.fasta.out",        
        trfMasked="trf/to_mask.{index}.fasta.trf"
    output:
        comb="comb/to_mask.{index}.fasta.masked",
        combOut="comb/to_mask.{index}.fasta.out",
    resources:
        mem = 8,		
        threads = 1,
        hrs = 4
    params:
        sd=SD
    shell:"""
export LD_LIBRARY_PATH=/home/hsiehph/shared/software/packages/htslib && {params.sd}/comask {output.comb} {input.humLib} {input.t2tLib} {input.trfMasked}
echo {input.humLibOut} > comb/to_mask.{wildcards.index}.names
echo {input.t2tLibOut} >> comb/to_mask.{wildcards.index}.names
{params.sd}/AppendOutFile.py {output.combOut} {params.sd}/repeat_masker.out.header comb/to_mask.{wildcards.index}.names
"""

rule CombineMaskedFasta:
    input:
        maskedFasta=lambda wildcards:  expand("comb/to_mask.{{name}}_{coords}.fasta.masked", coords=splitByContig[wildcards.name])
    output:
        combMaskedFasta="comb/masked.{name}.fasta"
    resources:
        mem = 16,
        threads = 1,
        hrs = 24
    params:
        sd=SD,
    shell:"""
module load python3/3.9.3_anaconda2021.11_mamba && {params.sd}/CombineFasta.py {output.combMaskedFasta} {wildcards.name} {SplitOverlap} {input.maskedFasta}
"""    
        
rule WriteOutNames:
    input:
        maskedContigsOut=expand("comb/to_mask.{index}.fasta.out", index=splitRegions),
    output:
        contigOutFOFN="comb/comb.out.fofn",
    resources:
        mem = 8,		
        threads = 1,
        hrs = 4
    run:
        outFile=open(output.contigOutFOFN,'w')
        for i in input.maskedContigsOut :
            outFile.write(i + "\n")
        outFile.close()



rule CatMaskedFasta:
    input:
        maskedContigs=expand("comb/masked.{name}.fasta", name=contigNames)
    output:
        maskedGenome=szMaskedAssembly,
    resources:
        mem = 8,		
        threads = 1,
        hrs = 4
    shell:"""
cat {input.maskedContigs} > {output.maskedGenome}
"""

rule bgzipRepeatMaskedAssembly:
	input: szMaskedAssembly
	output: szCompressedMaskedAssembly
	resources:
		mem = 8,		
		threads = 1,
		hrs = 4
	shell: """
		module load htslib/1.17-19-g07638e1 && bgzip {szMaskedAssembly}
"""

rule CombineMask:
    input:
        contigOutFOFN="comb/comb.out.fofn"
    output:
        maskedGenomeOut=szLocationsOfRepeats
    resources:
        mem = 8,		
        threads = 1,
        hrs = 4
    params:
        sd=SD
    shell:"""
{params.sd}/AppendOutFile.py {output.maskedGenomeOut} {params.sd}/repeat_masker.out.header  {input.contigOutFOFN}
"""
