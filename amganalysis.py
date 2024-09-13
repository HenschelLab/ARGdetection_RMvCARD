from math import pi
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [14, 20]
import pandas as pd
pd.set_option("display.max_columns", None)
from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

def loadRGI(sid, mapqThr = 10, pCovThr = 90, mappedReadThr = 10):
    gene = pd.read_csv(sid, sep='\t')
    gene['Reference Length'] = [int(str(rl).split(';')[0]) for rl in gene['Reference Length']]
    gene['coverage'] = 150*gene['Completely Mapped Reads'].div(gene['Reference Length'])  
    gene["relAb"] = 100*gene.coverage/gene.coverage.sum()
    # QC
    geneQC = gene[(gene['Average MAPQ (Completely Mapped Reads)'] > mapqThr) & (gene['Average Percent Coverage'] > pCovThr) & (gene['Completely Mapped Reads'] > mappedReadThr)]
    return geneQC

def extractHeatmapData(sid, col1='ARO Term', col2='Completely Mapped Reads', mapqThr = 10, pCovThr = 95, mappedReadThr = 10, log=False):
    geneQC = loadDF(sid, mapqThr=mapqThr, pCovThr=pCovThr, mappedReadThr=mappedReadThr)
    if log:
        geneQC[col2] = np.log(geneQC[col2])
    return dict(zip(geneQC[col1], geneQC[col2]))

def cluster(sampleFiles, groupbyCol='AMR Gene Family', **kwargs):
    '''Clustering analysis for RGI output files
    This function assumes that the user provides a list of files produced by a run of RGI (commonly ending in 'gene_mapping_data.txt').
    The clustering can happen by any appropriate column name from those files.
    Basic QC is applied using the loadRGI function.
    '''
    sampleIDs = [sid.split('/')[-1].split('.')[0] for sid in sampleFiles]
    card = pd.concat([loadRGI(sid, **kwargs).groupby(groupbyCol).count().iloc[:,:1] for sid in sampleFiles], axis=1)
    card.columns = sampleIDs #[f'{loctype[:3]}/{citycode[loc]}/AFY{sid}' for loctype, loc, sid in zip(loctypes, locs, sampleIDs)]
    clustermap = sns.clustermap(card.corr())
    
    ## heatmap
    orderedSamples = [label.get_text() for label in clustermap.ax_heatmap.get_xticklabels()]
    fig, axs = plt.subplots (figsize=(10, 21))
    sns.heatmap(card[orderedSamples], annot=True)