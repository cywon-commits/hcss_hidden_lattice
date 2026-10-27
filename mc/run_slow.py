"""Slow NVT cooling at rho*=0.98 (HCSS). Finer T steps near ordering, more sweeps per T."""
import sys, os, time, math, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "engine"))
import hcss_mc as mc
lam=float(sys.argv[1]); out=sys.argv[2]; nsw=int(sys.argv[3]); seed=int(sys.argv[4])
nx,ny,rho=32,36,0.98
a=math.sqrt(2/(math.sqrt(3)*rho))
pts=np.array([[(i+0.5*(j%2))*a, j*a*math.sqrt(3)/2] for j in range(ny) for i in range(nx)])
Lx,Ly=nx*a, ny*a*math.sqrt(3)/2; box=np.array([Lx,0.0,Ly]); s=np.stack([pts[:,0]/Lx,pts[:,1]/Ly],1)%1.0
n=mc.total_count(s,*box,lam); assert n>=0
temps=[0.6,0.45]+[round(x,3) for x in np.arange(0.36,0.195,-0.02)]
dmax=0.05; t0=time.time()
for T in temps:
    beta=1/T; acc=0; es=[]
    for k in range(nsw):
        n,ap,ab=mc.sweep(s,box,lam,beta,0.0,dmax,0.0,0.0,n,seed); seed+=1; acc+=ap
        if (k+1)%2000==0:
            r=acc/(2000*len(s)); acc=0; dmax*=1.1 if r>0.4 else 0.9; dmax=min(max(dmax,0.005),0.3)
        if k>=nsw//2 and k%500==0: es.append(n/len(s))
    print(json.dumps({"T":T,"E":float(np.mean(es)),"t":round(time.time()-t0)}),flush=True)
    np.savez(f"{out}_T{T:.3f}.npz",s=s,box=box,lam=lam,T=T)
