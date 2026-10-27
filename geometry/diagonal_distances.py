import numpy as np
from scipy.spatial import cKDTree
from sqtri import *
src=open('develop.py').read().split('if __name__')[0].replace("return maxerr,dmin,len(Q)","return Q,maxerr")
ns={}; exec(src,ns)
rng=np.random.default_rng(3)
P=DOD[0].copy()
for k in range(2): P=inflate(P,lambda i,p:int(rng.integers(2)))
nbr,tris,sqs=build(P); I=interior_vertices(P,nbr,tris,sqs)
sigma,*_=cover(P,nbr,sqs,I)
R=np.linalg.norm(P,axis=1).max()
keep=lambda f: all(np.linalg.norm(P[v])<0.75*R for v in f)
ti=[f for f in tris if keep(f)]; idx=[s for s,f in enumerate(sqs) if keep(f)]
si=[sqs[s] for s in idx]; sg={k:sigma[s] for k,s in enumerate(idx)}
core=[v for v in set(v for f in ti+si for v in f) if np.linalg.norm(P[v])<0.45*R]   # count pairs from deep vertices only
Nt,Ns=len(tris),len(sqs); N=Nt/2+Ns
def area_pp(th):  # per particle, a=1
    return (Nt*np.sqrt(3)/4+Ns*np.cos(np.radians(th)))/N
thetas=np.arange(0,60.01,0.5)
D={}
for th in thetas:
    Q,err=ns['develop'](P,ti,si,sg,th)
    vs=sorted(Q); X=np.array([Q[v] for v in vs]); pos={v:i for i,v in enumerate(vs)}
    t=cKDTree(X); ci=[pos[v] for v in core]
    d,_=t.query(X[ci],k=40)
    D[th]=d[:,1:]   # distances (a=1 units) from deep vertices
np.savez('vtheta_dist.npz',thetas=thetas,D=np.array([D[t] for t in thetas]),Nt=Nt,Ns=Ns)
print("saved; closure ok; deep vertices",len(core))
