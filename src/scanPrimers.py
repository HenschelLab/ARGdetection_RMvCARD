#scanPrimers.py
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import pylcs
from collections import Counter
card = "card.fasta"

primers = pd.read_excel('../Resistome/primers2.0.xlsx')

def searchPrimerPair(forward, reverse, card):
    for rec in SeqIO.parse(card, 'fasta'):
        reverseComp = Seq(reverse).reverse_complement()
        forwardComp = Seq(forward).reverse_complement()
        if forward in rec.seq and reverseComp in rec.seq:
            yield rec.id
        elif forwardComp in rec.seq and reverse in rec.seq:
            yield rec.id + " [RVS]"

def shortID(matches, field=2): return [m.split('|')[field].split(':')[-1] for m in matches]
def getCommon(matches, thresh=100):
    '''identify the part that all genes have in common'''
    if not matches: return ''
    gdf = pd.DataFrame([list(m) for m in matches])
    common = ''
    for col in gdf.columns.values:
        c = Counter(gdf[col])        
        if len(c) == 1: common += list(c.keys())[0]
        else: break
    return common

def lcstext(a,b):
    res = pylcs.lcs_sequence_idx(a, b)
    return ''.join([b[i] for i in res if i!= -1])

matches2rep = 500
header = 'Suggestion_Gene target_AMG Class_Nr of Matches'.split('_')
header += ['match %i'%i for i in range(1,11)] 
data = []
nameMatchTable = []
nameMatchTableHeader = ['RM', 'CARDcons', 'lcs(Rm,Card)', 'lcs(rm,card)'] + ['match %i'%i for i in range(1,matches2rep + 1)]

for idx, row in primers.iterrows():
    matches = [match for match in searchPrimerPair(row["Forward primer"], row["Reverse Primer"], card)]
    matches1 = shortID(matches)
    matchesARO = shortID(matches, 0)

    nameCons100 = getCommon(matches1)
    #nameCons90 = getCommon(matches1, thresh=90) ## not yet implemented
    lcs = lcstext(row['Gene target'], nameCons100)
    lcs2 = lcstext(row['Gene target'].lower(), nameCons100.lower())
    data.append([row["Suggestion "], row["Gene target"], row['Unnamed: 4'], len(matches) ] + matches[:10] + [nameCons100, lcs])

    if len(matches)>-1:
        nameMatchTable.append([row['Gene target'], nameCons100, lcs, lcs2] + matchesARO[:matches2rep])
#    if len(data) > 50:
#        break
ndf = pd.DataFrame(nameMatchTable) 
ndf.columns=nameMatchTableHeader[:ndf.shape[1]]
ndf.to_csv('primerCrossmatch2all.csv') 
#pd.DataFrame(data, columns=header).to_excel('primerCrossmatch2.xlsx')
