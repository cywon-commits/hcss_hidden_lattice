import numpy as np
from collections import Counter
from scipy.spatial import cKDTree
from sqtri import *
s3=np.sqrt(3); c=(1+s3)/np.sqrt(2)

def lattice(basis, motif, R):
    a1,a2=map(np.array,basis); pts=[]
    n=int(R/min(np.linalg.norm(a1),np.linalg.norm(a2)))+3
    for i in range(-n,n+1):
        for j in range(-n,n+1):
            for m in motif:
                p=i*a1+j*a2+np.array(m)
                if np.linalg.norm(p)<R: pts.append(p)
    return dedupe(np.array(pts))

def sq_corners(center,rot):
    return [center+np.array([np.cos(a),np.sin(a)])/np.sqrt(2) for a in np.radians(rot+45+90*np.arange(4))]
seeds={
 "tri 3^6": ([(1,0),(0.5,s3/2)],[(0,0)]),
 "sq 4^4": ([(1,0),(0,1)],[(0,0)]),
 "elong 3^3.4^2": ([(1,0),(0.5,1+s3/2)],[(0,0),(0,1)]),
 "snub 3^2.4.3.4": ([(c,0),(0,c)], [tuple(p) for p in sq_corners(np.zeros(2),15)+sq_corners(np.array([c/2,c/2]),-15)]),
}

def holonomy(name,P,periods):
    nbr,tris,sqs=build(P)
    I=interior_vertices(P,nbr,tris,sqs)
    sigma,conf,comps,edges,vsq=cover(P,nbr,sqs,I)
    cen=np.array([P[f].mean(0) for f in sqs]); tree=cKDTree(cen)
    r=np.linalg.norm(P,axis=1).max()
    deg=Counter(len(nbr[v]) for v in I)
    out=[]
    for T in periods:
        T=np.array(T); same=flip=0
        for s,f in enumerate(sqs):
            if np.linalg.norm(cen[s])>0.4*r: continue
            d,t=tree.query(cen[s]+T)
            if d>1e-4: continue
            a,b=(f[0],f[2]) if sigma[s]==1 else (f[1],f[3])
            g=sqs[t]; a2,b2=(g[0],g[2]) if sigma[t]==1 else (g[1],g[3])
            pa={tuple(np.round(P[a]+T,4)),tuple(np.round(P[b]+T,4))}
            pb={tuple(np.round(P[a2],4)),tuple(np.round(P[b2],4))}
            if pa==pb: same+=1
            else: flip+=1
        out.append((tuple(np.round(T,3)),same,flip))
    print(f"[{name}] N_int={len(I)} deg={dict(sorted(deg.items()))} sq={len(sqs)} Qcomps={comps} conflicts={len(conf)}")
    for T,sm,fl in out: print(f"   T={T}: same={sm} flip={fl}  -> holonomy {'+1' if fl==0 else ('-1' if sm==0 else 'MIXED')}")

for name,(B,M) in seeds.items():
    if name not in ("tri 3^6","sq 4^4"):
        P=lattice(B,M,30); holonomy(name,P,B)
    for orient in ("uniform","alt"):
        P0=lattice(B,M,14)
        # alt: orientation alternates on a 2x supercell (period 2B)
        if orient=="uniform": ch=lambda i,p:0; per=[LAM*np.array(b) for b in B]
        else:
            A=np.array(B).T; Ainv=np.linalg.inv(A)
            def ch(i,p):
                k=np.round(Ainv@p,6); return int((np.floor(k[0]+1e-6)+2*np.floor(k[1]+1e-6))%2) if len(M)==1 else int(np.floor(k[0]+1e-6)%2)
            per=[LAM*np.array(B[0])*2, LAM*np.array(B[1])*2]
        P=inflate(P0,ch); holonomy(f"infl({name}) {orient}",P,per)
