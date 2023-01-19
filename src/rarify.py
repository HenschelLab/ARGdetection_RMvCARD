## conda environment bio
"""
In order to quantify, how comprehensive the detection of AMGs, but also drug classes and 
Rarefaction was conducted by randomly subsampling reads from the complete sample AFY01. 
This procedure yields forward and reverse fastq files for subsets of 10%, 20% ... 90% of reads, with 3 replica respectively. 

"""
from Bio import SeqIO
import numpy as np
import gzip

def newFilename(filename, fraction, rep=0):
    fn = filename.split('.')
    ext = '.'.join(fn[1:])
    return f"{fn[0]}_{fraction:02d}_{rep:02}.{ext}"

## to be caught with glob
filename = 'Data/AFY01_1.fastq.gz'
filename2 = 'Data/AFY01_2.fastq.gz'

i1 = SeqIO.parse(gzip.open(filename,'rt'), 'fastq')
i2 = SeqIO.parse(gzip.open(filename2,'rt'), 'fastq')

out1s = {}
out2s = {}
## creating an array of open file handles
for fraction0 in range(1,10):
    fraction = fraction0/10
    for rep in range(3):
        newFn1 = newFilename(filename, fraction0, rep)
        newFn2 = newFilename(filename2, fraction0, rep)
        out1s[(fraction0, rep)] = gzip.open(newFn1,'wt')
        out2s[(fraction0, rep)] = gzip.open(newFn2,'wt')

## The actual rarefaction
## Datasets are created simultaneously
## To keep the memory footprint low, we immediately write to disk, based on probability.
## Also, of note, we subsample forward and reverse reads in pairs.
for rec1, rec2 in zip(i1, i2):
    for rep in range(3):
        for fraction0 in range(1,10):
            fraction = fraction0/10
            if np.random.rand() < fraction: 
                out1 = out1s[(fraction0, rep)]
                out2 = out2s[(fraction0, rep)]
                SeqIO.write(rec1, out1, 'fastq')
                SeqIO.write(rec2, out2, 'fastq')

## Closing all open files
for out1 in out1s.values():
    out1.close()
for out2 in out2s.values():
    out2.close()

