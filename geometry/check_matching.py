import numpy as np, networkx as nx
from sqtri import *
rng=np.random.default_rng(2)
P=DOD[0].copy()
for k in range(2): P=inflate(P, lambda i,p:int(rng.integers(2)))
nbr,tris,sqs=build(P); I=interior_vertices(P,nbr,tris,sqs)
# deep interior: vertices whose all neighbours are interior
D={v for v in I if all(u in I for u in nbr[v])}
odd=[v for v in nbr if len(nbr[v])==5 and v in I]
G=nx.Graph(); G.add_nodes_from(odd)
for v in odd:
    for u in nbr[v]:
        if u in G: G.add_edge(u,v)
M=nx.max_weight_matching(G,maxcardinality=True)
matched=set(); [matched.update(e) for e in M]
um=[v for v in odd if v not in matched]
umd=[v for v in um if v in D and all(u in D for u in nbr[v])]
r=np.linalg.norm(P,axis=1).max()
print("verts",len(P),"odd interior",len(odd),"matched",len(matched),"unmatched",len(um),
      "unmatched deep-interior",len(umd))
print("radii of unmatched / patch radius:",np.round(sorted(np.linalg.norm(P[um],axis=1)/r)[:10],2))
