"""NPT compression series for HCSS (delta=1.27) at fixed T: sweeps per pressure, snapshot after each step."""
import sys, os, json, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "engine"))
import hcss_mc as mc
src,T,out,nsw,seed=sys.argv[1],float(sys.argv[2]),sys.argv[3],int(sys.argv[4]),int(sys.argv[5])
Ps=[float(x) for x in sys.argv[6].split(",")]
z=np.load(src); s=z["s"].copy(); box=z["box"].copy(); lam=float(z["lam"])
n=mc.total_count(s,*box,lam); assert n>=0
dmax,dbox,dsh=0.03,1e-4,1e-3; t0=time.time(); N=len(s)
for P in Ps:
    acc=accb=0; rho=[]; e=[]
    for k in range(1,nsw+1):
        n,ap,ab=mc.sweep(s,box,lam,1/T,P,dmax,dbox,dsh,n,seed); seed+=1; acc+=ap; accb+=ab
        if k%200==0:
            r=acc/(200*N); rb=accb/200; acc=accb=0
            dmax*=1.1 if r>0.4 else 0.9; dmax=min(max(dmax,0.003),0.3)
            dbox*=1.1 if rb>0.35 else 0.9; dbox=min(max(dbox,1e-6),0.02); dsh=dbox*10
        if k>nsw//2 and k%200==0: rho.append(N/(box[0]*box[2])); e.append(n/N)
    np.savez(f"{out}_P{P:g}.npz",s=s,box=box,lam=lam,T=T,P=P)
    print(json.dumps({"P":P,"rho":round(float(np.mean(rho)),4),"rho_sd":round(float(np.std(rho)),4),"E":round(float(np.mean(e)),4),"t":round(time.time()-t0)}),flush=True)
