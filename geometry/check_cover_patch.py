import numpy as np
from collections import Counter
from sqtri import *
rng = np.random.default_rng(1)
def test(name, P):
    nbr, tris, sqs = build(P)
    I = interior_vertices(P, nbr, tris, sqs)
    sigma, conf, comps, edges, vsq = cover(P, nbr, sqs, I)
    deg = {v: len(nbr[v]) for v in I}
    N = len(I); n5 = sum(1 for v in I if deg[v]==5)
    Nt, Ns = len(tris), len(sqs)
    x = Nt/Ns if Ns else np.inf
    d2 = degrees_after(P, nbr, sqs, sigma)
    rs = {s: int(rng.choice([1,-1])) for s in range(len(sqs))}
    d3 = degrees_after(P, nbr, sqs, rs)
    print(f"[{name}] verts={len(P)} interior={N} tris={Nt} sq={Ns} x={x:.4f}")
    print(f"  interior degrees {dict(sorted(Counter(deg.values()).items()))}")
    print(f"  n_odd/N={n5/N:.4f}  formula={2/(1+x/2):.4f}  (12-fold 4√3-6={4*np.sqrt(3)-6:.4f})")
    print(f"  Q comps={comps}  conflicts={len(conf)}")
    print(f"  after cover: interior deg dist {dict(sorted(Counter(d2[v] for v in I).items()))}")
    print(f"  random diagonals n_odd/N={sum(d3[v]%2 for v in I)/N:.4f} pred={n5/(2*N):.4f}")
P = DOD[0].copy()
Pd = P
for k in range(3): Pd = inflate(Pd, lambda i,p:0)
test("Stampfli all-0 x3", Pd)
Pr = P
for k in range(3): Pr = inflate(Pr, lambda i,p: int(rng.integers(2)))
test("random orient x3", Pr)
