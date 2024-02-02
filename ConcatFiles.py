#!/usr/bin/env python
import sys
inFofn=sys.argv[1]
outFilename=sys.argv[2]

outFile=open(outFilename)
for i in inFofn:
    f = open(i)
    outFile.write(f.readlines())
    f.close()
