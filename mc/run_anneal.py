import sys, os, json, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "engine"))
import hcss_mc as mc
src,T,out,nsw,every,seed=sys.argv[1],float(sys.argv[2]),sys.argv[3],int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6])
z=np.load(src); s=z["s"].copy(); box=z["box"].copy(); lam=float(z["lam"])
n=mc.total_count(s,*box,lam); dmax=0.03; acc=0; t0=time.time()
for k in range(1,nsw+1):
    n,ap,ab=mc.sweep(s,box,lam,1/T,0.0,dmax,0.0,0.0,n,seed); seed+=1; acc+=ap
    if k%2000==0:
        r=acc/(2000*len(s)); acc=0; dmax*=1.1 if r>0.4 else 0.9; dmax=min(max(dmax,0.005),0.3)
    if k%every==0:
        np.savez(f"{out}_{k//1000:05d}k.npz",s=s,box=box,lam=lam,T=T)
        print(json.dumps({"k":k,"E":n/len(s),"t":round(time.time()-t0)}),flush=True)
