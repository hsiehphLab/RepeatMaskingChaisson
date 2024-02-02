import os
import tempfile
import subprocess
import os.path


# Config
configfile: "t2t_repeat_mask.json"

if "asm" not in config and "assembly" in config:
    config["asm"] = config["assembly"]

assembly="assembly.orig.fasta"
asmFai="assembly.orig.fasta.fai"
#asmFai=config["asm"] + ".fai"
asmFaiFile=open(asmFai)
faiLines=asmFaiFile.readlines()

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
        mask="assembly.repeat_masked.fasta",
        maskedGenomeOut="assembly.repeat_masked.fasta.out"        


rule IndexGenome:
    input:
        asm="assembly.orig.fasta"
    output:
        fai="assembly.orig.fasta.fai"
    params:
        grid_opts=config["grid_small"]
    shell:"""
samtools faidx {input.asm}
"""

rule SplitGenome:
    input:
        asm="assembly.orig.fasta",
        fai="assembly.orig.fasta.fai",
    output:
        split=expand("split/to_mask.{region}.fasta", region=splitRegions)
    params:
        grid_opts=config["grid_medium"],
        sd=SD
    resources:
        load=1
    shell:"""
mkdir -p split
{params.sd}/DivideFasta.py {input.asm} {SplitSize} {SplitOverlap} split/to_mask
"""

rule MaskContig:
    input:
        contig="split/to_mask.{index}.fasta"
    output:
        mask="masked/to_mask.{index}.fasta.masked",
        maskOut="masked/to_mask.{index}.fasta.out"
    conda:
        "rmsk"
    params:
        grid_opts=config["grid_repeatmasker"],
        repeatLibrary=config["repeat_library"],
        sd=SD        
    shell:"""
mkdir -p masked
TEMP="$TMPDIR/$$_$RANDOM/"
mkdir -p $TEMP
cp \"{input.contig}\" \"$TEMP/to_mask.{wildcards.index}.fasta\" && \
pushd $TEMP &&  \
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
    params:
        grid_opts=config["grid_repeatmasker"],
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
  {params.sd}/hardmask {input.mask} $TEMP/to_mask.{wildcards.index}.fasta && \
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
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
mkdir -p trf    
{params.sd}/trf-mod -p 500 -l 20000 {input.orig} > {output.trf}.bed
{params.sd}/bemask {input.orig} {output.trf}.bed {output.trf}
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
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
{params.sd}/comask {output.comb} {input.humLib} {input.t2tLib} {input.trfMasked}
echo {input.humLibOut} > comb/to_mask.{wildcards.index}.names
echo {input.t2tLibOut} >> comb/to_mask.{wildcards.index}.names
{params.sd}/AppendOutFile.py {output.combOut} {params.sd}/repeat_masker.out.header comb/to_mask.{wildcards.index}.names
"""

rule CombineMaskedFasta:
    input:
        maskedFasta=lambda wildcards:  expand("comb/to_mask.{{name}}_{coords}.fasta.masked", coords=splitByContig[wildcards.name])
    output:
        combMaskedFasta="comb/masked.{name}.fasta"
    params:
        sd=SD,
        grid_opts=config["grid_medium"]
    shell:"""
{params.sd}/CombineFasta.py {output.combMaskedFasta} {wildcards.name} {SplitOverlap} {input.maskedFasta}
"""    
        
rule WriteOutNames:
    input:
        maskedContigsOut=expand("comb/to_mask.{index}.fasta.out", index=splitRegions),
    output:
        contigOutFOFN="comb/comb.out.fofn",
    params:
        grid_opts=config["grid_small"]
    run:
        outFile=open(output.contigOutFOFN,'w')
        for i in input.maskedContigsOut :
            outFile.write(i + "\n")
        outFile.close()



rule CatMaskedFasta:
    input:
        maskedContigs=expand("comb/masked.{name}.fasta", name=contigNames)
    output:
        maskedGenome="assembly.repeat_masked.fasta",
    params:
        grid_opts=config["grid_small"]
    shell:"""
cat {input.maskedContigs} > {output.maskedGenome}
"""
    
    
rule CombineMask:
    input:
        contigOutFOFN="comb/comb.out.fofn"
    output:
        maskedGenomeOut="assembly.repeat_masked.fasta.out"        
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
{params.sd}/AppendOutFile.py {output.maskedGenomeOut} {params.sd}/repeat_masker.out.header  {input.contigOutFOFN}
"""
        
