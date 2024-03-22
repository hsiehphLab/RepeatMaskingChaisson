#!/usr/bin/env python

import sys
import re

szInputFileName = sys.argv[1]
# input file will look like this:
# to_mask.haplotype1-0000031_2380000_3470964_20000_1.fasta.out

szPattern = "to_mask\..+_(?P<bias>\d+)_(\d+)_(\d+)_(\d+).fasta.out$"

match = re.search( szPattern, szInputFileName )
if ( match ):
    nBias = int( match.group("bias" ) )
    print( "will use bias: " + match.group("bias" ) )
else:
    exit( "didn't match " + szInputFileName )


with open( szInputFileName, "r" ) as fInputMasked, open( sys.argv[2], "w" ) as fOutputMasked:
    nLine = 0
    while True:
        szLine = fInputMasked.readline()
        if ( szLine == "" ):
            break

        nLine += 1

        # copy over header:
        #   SW   perc perc perc  query               position in query           matching                                 repeat               position in repeat
        #score   div. del. ins.  sequence            begin   end        (left)   repeat                                   class/family      begin   end    (left)   ID
        #  (blank line)
        

        if ( nLine <= 3 ):
            fOutputMasked.write( szLine )
            continue

        aWords = szLine.split()

        # looks like:
        #   561   11.5  1.2  0.6  haplotype1-0000031        1    7409 (1083555) C SAR                                      Satellite            (11)   7453      1    1 *

        #    0     1     2    3     4                       5     6

        aWords[5] = str( int( aWords[5] ) + nBias )
        aWords[6] = str( int( aWords[6] ) + nBias )

        szOutputLine = "  ".join( aWords )

        fOutputMasked.write( szOutputLine + "\n" )




