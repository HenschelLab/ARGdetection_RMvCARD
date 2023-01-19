from collections import Counter
import pandas as pd
from goatools.obo_parser import GODag
import networkx as nx
from itertools import combinations

def parents(node):
    if node.id == "ARO:1000001": return [] ## adjust if multiple roots
    result = []
    for parent in node.parents:
        result += parents(parent)
    result += [node]
    return result

def lineage(node):
    if node.id == "ARO:1000001": return [] ## adjust if multiple roots
    result = []
    for parent in node.parents:
        result += parents(parent)
    result += [node]
    return result


parentIDset = lambda parentSet: {parent.id for parent in parentSet}

def intersectionAll(sets):
    if not sets: return set()
    isec= sets[0]
    for s in sets[1:]:
        isec = isec.intersection(s)
    return isec
    
def unionAll(sets):
    if not sets: return set()
    uni = sets[0]
    for s in sets[1:]:
        uni = uni.union(s)
    return uni

def countStats(nodes):
    nodes = flatten(nodes)
    stats = Counter(nodes)
    return stats.most_common(1)[0][1]/len(nodes)

def flatten(x): return [item for sublist in x for item in sublist]

def commonParents(nodes):
    parents = [parentIDset(node.parents) for node in nodes]
    return intersectionAll(parents), countStats(parents)

def findParents(nodes):
    return unionAll([parentIDset(node.parents) for node in nodes])

def commonGrandParents(nodes):
    grandparents = [findParents(node.parents) for node in nodes]
    return intersectionAll(grandparents), countStats(grandparents)

def mrca(nodes1, nodes2):
    p1 = unionAll([node.parents for node in nodes1])
    p2 = unionAll([node.parents for node in nodes2])
    commonParent = p1.intersection(p2)
    if commonParent: return commonParent
    return mrca(p1, p2)

def addEdge(G, node):
    for child in node.children:
        edge = node.id, child.id
        if not G.has_edge(*edge):
            G.add_edge(*edge)
            addEdge(G, child)

def makeGraph(g):
    import networkx as nx
    roots = [n for n in g.values() if len(n.parents)==0]
    G = nx.DiGraph()
    for root in roots:
        addEdge(G, root)
    return G

def levelDiffrep(c):
    (n1, n2), ancestor = c
    return min(gi[n1].level, gi[n2].level) - gi[ancestor].level 
    
def maxDiversity(nodes):
    pairs = combinations(nodes, 2)
    ancestry = [ (levelDiffrep(c),c) for c in nx.all_pairs_lowest_common_ancestor(G, pairs=pairs) ]
    return max(ancestry)
        

g = GODag("aro.obo")
G = makeGraph(g)

gd = {node.name:node for node in g.values() }
gi = {node.id:node for node in g.values() }
df = pd.read_csv('ResistomeComparison/primerCrossmatch2all2.csv')
addcols = []
colnames = 'matches comParentS comParPct comGrandParentS comGrandParPct max_div n1 n2 mrca'.split()
#df1 = df[df.RM.isin({'ampC', 'blaACC-1'})]
for i,row0 in df.iterrows():
    row = ['ARO:%d'%match for match in row0.iloc[5:].dropna()]
    if len(row)>1:
        nodes = [gi[aro] for aro in row]
        cp, cpPct = commonParents(nodes)
        cg, cgPct = commonGrandParents(nodes)
        maxDiv, ((n1, n2), ancestor) = maxDiversity(row)
        infoName = gi[n1].name, gi[n2].name, gi[ancestor].name
        infoLevel = gi[n1].level, gi[n2].level, gi[ancestor].level
        addcols.append((len(row), cp, cpPct, cg, cgPct, maxDiv, n1, n2, ancestor, *infoName, *infoLevel))
    else:
        addcols.append( (0, None, None, None, None, None, None, None, None, None, None, None, None, None) )

parentStats = pd.DataFrame(addcols, columns=colnames)
primerCross = pd.concat([df.iloc[:,:5], parentStats], axis=1)

#primerCross[primerCross.matches>1].sort_values(by='comParPct').to_csv('ResistomeComparison/primerCrossmatch2all2_homologyStats.csv')
