#!/usr/bin/env python

import os
import shutil

nFilesToKeepInDirectories = 5


aDirectoriesToMainlyEmpty = ["masked", "comb", "split", "t2t", "trf", "comb_fixed_coords" ]


for szDirectory in aDirectoriesToMainlyEmpty:
    print( "now working on " + szDirectory )

    nFilesSoFar = 0

    with os.scandir( szDirectory) as it:
        for entry in it:
            if not entry.name.startswith('.') and entry.is_file():
                nFilesSoFar += 1
                if ( nFilesSoFar > nFilesToKeepInDirectories ):
                    #print( "about to delete: " + entry.name )
                    os.remove( szDirectory + "/" + entry.name )
                else:
                    print( "will keep: " + entry.name )




# delete logs.  This saves each subdirectory of logs and saves 5 logs
# from each subdirectory.
szDirectory = "logs"
with os.scandir( szDirectory ) as it:
    for entry in it:
        if not entry.name.startswith('.') and entry.is_dir():

            nFilesSoFar = 0

            szSubDir = szDirectory + "/" + entry.name
            with os.scandir( szSubDir ) as subdirEntry:
                for ssubdirEntry in subdirEntry:

                    if not ssubdirEntry.name.startswith('.') and ssubdirEntry.is_file():
                        nFilesSoFar += 1
                        if ( nFilesSoFar > nFilesToKeepInDirectories ):
                            #print( "about to delete " + szSubDir + "/" + ssubdirEntry.name )
                            os.remove( szSubDir + "/" + ssubdirEntry.name )
                        else:
                            print( "will keep: " + szSubDir + "/" + ssubdirEntry.name )


szDirectory = ".snakemake/metadata"
if ( os.path.isdir( szDirectory )):
     print( "will delete: " + szDirectory )
     shutil.rmtree( szDirectory )
else:
    print( szDirectory + " is already deleted" )

szFileToDelete = "masked_assembly_without_windowmasker.fa"    
if os.path.exists( szFileToDelete ):
    os.remove( szFileToDelete )    

shutil.rmtree( "windowmasker", ignore_errors = True )
