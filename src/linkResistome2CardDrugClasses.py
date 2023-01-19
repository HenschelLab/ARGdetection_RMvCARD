## conda activate goatools
from collections import Counter
import pandas as pd
from goatools.obo_parser import GODag
import re
import pdb

"""This script takes gene mapping files as produced by rgi bwt and derives a secondary annotation that uses terms that are used in ResistoMap.
The algorithm first establishes a top-level mapping between CARD identifiers and equivalent drug classes used in ResistoMap.
For each identified gene in the CARD report, the algorithm checks, if the ARO association is directly in the top-level mapping. If not,
the algorithm recursively checks, if any ontological superclass of the gene can be mapped to ResistoMap. This was done with the help of the OBO parser from the python library GOATOOLS [Klopfenstein et al], operating on the see `linkResistome2CardDrugClasses.py`.

Klopfenstein DV, Zhang L, Pedersen BS, ... Tang H GOATOOLS: A Python library for Gene Ontology analyses Scientific reports | (2018) 8:10872 | DOI:10.1038/s41598-018-28948-z

If a term in CARD annotation is either directly listed in this map or it's a child of a term, it
The output is a csv file (ngs2resi...) which is further processed in rgi_visual2.ipynb, including barplots of antibiotic  
"""

g = GODag("aro.obo")
gd = {node.name:node for node in g.values() } ## lookup with proper name like "fluoroquinolone antibiotics"
root = [n for (nid, n) in g.items() if n.level==0][0]
#mapping terms from the Antibiotic resistance Ontology to terms that are used in ResistoMap:
#Rationale: if a term in CARD annotation is either directly listed in this map or it's a child of a term, we can map
resistome={"ARO:3000282": "sulfonamide antibiotic",
           "ARO:0000016": "aminoglycoside antibiotic",
           "ARO:3000387": "phenicol antibiotic",
           "ARO:3000081": "glycopeptide antibiotic", ## vancomycin
           "ARO:0000001": "fluoroquinolone antibiotic",
           "ARO:3000007": "beta-lactam antibiotic",
           "ARO:3000050": "tetracycline antibiotic",
           "ARO:3000188": "trimethoprim", ## or 3000171 (parent: diaminopyrimidine)
           "ARO:3000171": "diaminopyrimidine (trimethoprim)", ## or 3000171 (parent: diaminopyrimidine)
           "ARO:0000000": "mlsb", #"macrolide antibiotic",
           "ARO:0000017": "mlsb", # "lincosamide antibiotic",
           "ARO:0000026": "mlsb"} # streptogramin antibiotic"

resistomeIDs = set(resistome.keys())

def parents(node):
    if node.id == "ARO:1000001": return [] ## adjust if multiple roots
    result = []
    for parent in node.parents:
        result += parents(parent)
    result += [node]
    return result

def readGeneTable(sid):
    gene = pd.read_csv(f'ResistomeComparison/NoWildAFY{sid}.gene_mapping_data.txt', sep='\t')
    cols = list(gene.columns.values)
    gene['Reference Length'] = [int(rl.split(';')[0]) for rl in gene['Reference Length']]
    gene['coverage'] = 150*gene['Completely Mapped Reads'].div(gene['Reference Length'])
    ## QC
#    gene = gene[(gene.coverage>1) & (gene['Completely Mapped Reads'] > 10)]
    gene = gene[gene['Average Percent Coverage'] > 95]
    return gene


def getTally(sid):
    gene = readGeneTable(sid)
    misc = Counter()
    tally = Counter()
    for i,row in gene.iterrows():
        drugclass = row['Drug Class']
        coverage = row['coverage']
        equiResistome = []
        drugclasses = [drug.strip() for drug in drugclass.split(';') if drug.strip()]
        for drug in drugclasses:
            drug = gd[drug]
            parentIds = [p.id for p in parents(drug)]
            linkIDs = set(parentIds).intersection(resistomeIDs)
            links = set([resistome[linkID] for linkID in linkIDs])
            if len(links) == 0: ## no link to resistome categories
                misc += Counter(parentIds)        
            elif len(links) == 1:
                equiResistome.append(list(links)[0])
        value2add = 1 #coverage
        if not equiResistome:           ## unlinked, register as "others"
            tally['Others'] += value2add
        elif len(set(equiResistome)) > 1:
            tally['MDR'] += value2add           
        else: # unique mapping
            for resi in set(equiResistome): ## deals with unique and ambiguously mapped drug classes
                tally[resi] += value2add
    return tally

sampleIDs = '01 03 04 07 10 11 13 14 15 16 17 18 19 20 22 99'.split()
result = pd.DataFrame([getTally(sid) for sid in sampleIDs])
result.index = sampleIDs
#result.to_csv("ResistomeComparison/ngs2resiQC2.csv")

"""
Most common antibiotic drug classes not linked to Resistome classes
ARO:3000053	level-02	depth-02	peptide antibiotic [antibiotic_resistance]                        Count:  69
ARO:3005386	level-02	depth-02	disinfecting agents and antiseptics [antibiotic_resistance]       Count:  65
ARO:3000157	level-02	depth-02	rifamycin antibiotic [antibiotic_resistance]                      Count:  48
ARO:3000003	level-02	depth-02	antibiotic without defined classification [antibiotic_resistance] Count:  46
ARO:3000103	level-02	depth-02	aminocoumarin antibiotic [antibiotic_resistance]                  Count:  30
ARO:0000025	level-03	depth-03	fosfomycin [antibiotic_resistance]                                Count:  15
ARO:3000670	level-02	depth-02	pleuromutilin antibiotic [antibiotic_resistance]                  Count:  11
ARO:3000520	level-03	depth-03	isoniazid [antibiotic_resistance]                                 Count:  10

##MLSb (deprecated, now taken care of by 
           "ARO:0000000": "mlsb", #"macrolide antibiotic",                                             
           "ARO:0000017": "mlsb", # "lincosamide antibiotic",                                          
           "ARO:0000026": "mlsb"} # streptogramin antibiotic" 
id: ARO:3000347 name: ErmA
id: ARO:3000375 name: ErmB
id: ARO:3000392 name: Erm(37)
id: ARO:3000495 name: ErmD
id: ARO:3000498 name: ErmF
id: ARO:3000560 name: Erm 23S ribosomal RNA methyltransferase
id: ARO:3000593 name: ErmQ
id: ARO:3000595 name: ErmT
id: ARO:3000598 name: Erm(31)
id: ARO:3000599 name: Erm(33)
id: ARO:3000600 name: Erm(34)
id: ARO:3000601 name: Erm(38)
id: ARO:3000602 name: Erm(39)
id: ARO:3000603 name: Erm(41)
id: ARO:3000604 name: Erm(35)
id: ARO:3001265 name: Erm(30)
id: ARO:3003106 name: Erm(42)
id: ARO:3003908 name: Erm(47)
id: ARO:3005099 name: 23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(A)


MISC
    if len(set(equiResistome)) == 0:
        unlinked += 1
    if len(set(equiResistome)) == 1:
        unique += 1
    if len(set(equiResistome)) >1:
        ambiguous += 1
        print (drugclass, "-->", set(equiResistome))
        #pdb.set_trace()


"""
