import numpy as np
from collections import defaultdict, deque
from scipy.spatial import cKDTree
from sqtri import *

def develop(P, tris, sqs, sigma, theta_deg, return_positions=False):
    th=np.radians(theta_deg)
    # tile corner angles: tri 60; square corner covered by cover diagonal -> 90+th, else 90-th
    tiles=[]
    for f in tris: tiles.append((list(f),[np.pi/3]*3))
    for s,f in enumerate(sqs):
        cov=(f[0],f[2]) if sigma[s]==1 else (f[1],f[3])
        tiles.append((list(f),[np.pi/2+th if v in cov else np.pi/2-th for v in f]))
    # orient tiles ccw consistently with original positions
    for k,(f,a) in enumerate(tiles):
        x=P[f]; area=0.5*np.sum(x[:,0]*np.roll(x[:,1],-1)-np.roll(x[:,0],-1)*x[:,1])
        if area<0: tiles[k]=(f[::-1],a[::-1])
    edge_tiles=defaultdict(list)
    for k,(f,a) in enumerate(tiles):
        for i in range(len(f)): edge_tiles[frozenset((f[i],f[(i+1)%len(f)]))].append(k)
    Q={}; placed=set()
    f,a=tiles[0]
    Q[f[0]]=np.zeros(2); Q[f[1]]=np.array([1.0,0])
    def place(k):
        f,a=tiles[k]; n=len(f)
        # find a consecutive pair already placed
        for i in range(n):
            u,w=f[i],f[(i+1)%n]
            if u in Q and w in Q:
                d=Q[w]-Q[u]; ang=np.arctan2(d[1],d[0]); pos=Q[w]
                worst=0
                for j in range(1,n-1):
                    idx=(i+1+j-1)%n; nxt=f[(idx+1)%n]
                    ang+=np.pi-a[idx]; pos=pos+np.array([np.cos(ang),np.sin(ang)])
                    if nxt in Q: worst=max(worst,np.linalg.norm(Q[nxt]-pos))
                    else: Q[nxt]=pos
                return worst
        return None
    dq=deque([0]); placed={0}; maxerr=0
    while dq:
        k=dq.popleft(); e=place(k)
        if e is None: dq.append(k); continue
        maxerr=max(maxerr,e)
        f,_=tiles[k]
        for i in range(len(f)):
            for k2 in edge_tiles[frozenset((f[i],f[(i+1)%len(f)]))]:
                if k2 not in placed: placed.add(k2); dq.append(k2)
    if return_positions:
        return Q, maxerr
    X=np.array([Q[v] for v in sorted(Q)])
    dmin=cKDTree(X).query(X,k=2)[0][:,1].min()
    area=sum(np.sqrt(3)/4 for _ in tris)+sum(np.cos(np.radians(theta_deg)) for _ in sqs)
    return maxerr,dmin,len(Q)

if __name__=="__main__":
    rng=np.random.default_rng(3)
    P=DOD[0].copy()
    for k in range(2): P=inflate(P,lambda i,p:int(rng.integers(2)))
    nbr,tris,sqs=build(P); I=interior_vertices(P,nbr,tris,sqs)
    sigma,*_=cover(P,nbr,sqs,I)
    for th in (0,15,30,45,60,75):
        err,dmin,nv=develop(P,tris,sqs,sigma,th)
        print(f"theta={th:2d}: closure error={err:.2e}, min vertex distance={dmin:.3f}, vertices={nv}")
    # wrong cover (random diagonals) should fail to close
    rs={s:int(rng.choice([1,-1])) for s in range(len(sqs))}
    err,dmin,nv=develop(P,tris,sqs,rs,15)
    print(f"random diagonals theta=15: closure error={err:.2e}")
