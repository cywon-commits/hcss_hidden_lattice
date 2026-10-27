import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from collections import Counter
from sqtri import *
rng=np.random.default_rng(3)

def analyze(name,P,kdiag=0.0,eps=1e-6):
    nbr,tris,sqs=build(P); I=interior_vertices(P,nbr,tris,sqs)
    sigma,conf,comps,edges,vsq=cover(P,nbr,sqs,I)
    n=len(P); rows=[];cols=[];vals=[]
    def spring(i,j,k):
        e=P[j]-P[i]; e=e/np.linalg.norm(e); B=k*np.outer(e,e)
        for a,b,s in ((i,i,1),(j,j,1),(i,j,-1),(j,i,-1)):
            for x in range(2):
                for y in range(2):
                    rows.append(2*a+x);cols.append(2*b+y);vals.append(s*B[x,y])
    for i in nbr:
        for j in nbr[i]:
            if i<j: spring(i,j,1.0)
    if kdiag:
        for f in sqs: spring(f[0],f[2],kdiag); spring(f[1],f[3],kdiag)
    K=sp.csc_matrix((vals,(rows,cols)),shape=(2*n,2*n))+eps*sp.identity(2*n)
    lu=sla.splu(K)
    # eta_s = d|13| - d|02|  (>0 => diagonal 0-2 shorter => sigma=+1)
    def gvec(f):
        g=np.zeros(2*n)
        for (a,b),sg in (((f[1],f[3]),1),((f[0],f[2]),-1)):
            e=(P[b]-P[a]); e/=np.linalg.norm(e)
            g[2*b:2*b+2]+=sg*e; g[2*a:2*a+2]-=sg*e
        return g
    r=np.linalg.norm(P,axis=1).max(); cen=np.array([P[f].mean(0) for f in sqs])
    deep=[s for s in range(len(sqs)) if np.linalg.norm(cen[s])<0.5*r]
    deepset=set(deep)
    G={s:gvec(sqs[s]) for s in range(len(sqs))}
    X={s:lu.solve(G[s]) for s in deep}
    var={s:G[s]@X[s] for s in deep}
    agree=dis=0; rhos=Counter(); rlist=[]
    for s,t,J,v in edges:
        if s in deepset and t in deepset:
            rho=G[t]@X[s]/np.sqrt(var[s]*var[t])
            typ=len(nbr[v]); shared = sum(1 for u in sqs[s] if u in sqs[t])
            key="3^3.4^2 (shared edge)" if shared==2 else "3^2.4.3.4"
            rlist.append((key,rho,J))
    out={}
    for key,rho,J in rlist:
        d=out.setdefault(key,[0,0,[]]); d[0 if np.sign(rho)==J else 1]+=1; d[2].append(rho*J)
    print(f"[{name}] kdiag={kdiag}")
    for key,(a,b,l) in out.items():
        l=np.array(l); print(f"   {key}: sign(rho)=J in {a}, opposite in {b};  rho*J mean={l.mean():.3f} min={l.min():.3f} max={l.max():.3f}")
    # 2nd-neighbour squares (not Q-adjacent) correlation in the gauge: tau
    return sqs,cen,deep,X,G,var,sigma

P=DOD[0].copy()
for k in range(2): P=inflate(P,lambda i,p:int(rng.integers(2)))
for kd in (0.0,0.3):
    analyze("12-fold Stampfli x2",P,kd)
B=[(1,0),(0.5,np.sqrt(3)/2)]
