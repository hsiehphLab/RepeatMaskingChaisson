#!/usr/bin/env python


from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import sys
outFile=open(sys.argv[1], 'w')
contigName=sys.argv[2]
overlap=int(sys.argv[3])
inFiles = sys.argv[4:]

seq=""
first=True
for fn in inFiles:
    f = open(fn)
    for seqRec in SeqIO.parse(f, "fasta"):
        if first:
            seq+= str(seqRec.seq)
        else:
            seq+= str(seqRec.seq[overlap:])
        first=False
rec=SeqRecord.SeqRecord(Seq.Seq(seq), id=contigName, name="", description="")
SeqIO.write(rec, outFile, "fasta")
