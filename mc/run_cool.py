"""NVT cooling of a hexagonal lattice in the HCSS model (protocol of Pattabhiraman & Dijkstra 2017):
rho* = 0.98, cool T* step-wise at fixed density. Saves snapshots (fractional coords + box)."""
import sys, os, time, math, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "engine"))
import hcss_mc as mc

lam = float(sys.argv[1]); Tend = float(sys.argv[2]); out = sys.argv[3]
nsw = int(sys.argv[4]) if len(sys.argv) > 4 else 40000
nx, ny, rho = 32, 36, 0.98
a = math.sqrt(2 / (math.sqrt(3) * rho))
pts = np.array([[(i + 0.5 * (j % 2)) * a, j * a * math.sqrt(3) / 2] for j in range(ny) for i in range(nx)])
Lx, Ly = nx * a, ny * a * math.sqrt(3) / 2
box = np.array([Lx, 0.0, Ly])
s = np.stack([pts[:, 0] / Lx, pts[:, 1] / Ly], 1) % 1.0
n = mc.total_count(s, *box, lam)
assert n >= 0
temps = [t for t in (0.8, 0.6, 0.5, 0.45, 0.40, 0.36, 0.33, 0.30, 0.278, 0.26, 0.24, 0.22, 0.20) if t >= Tend - 1e-9]
dmax, seed = 0.05, 1
log = []
t0 = time.time()
for T in temps:
    beta = 1 / T
    acc = 0
    for k in range(nsw):
        n, ap, ab = mc.sweep(s, box, lam, beta, 0.0, dmax, 0.0, 0.0, n, seed); seed += 1
        acc += ap
        if (k + 1) % 2000 == 0:
            r = acc / (2000 * len(s)); acc = 0
            dmax *= 1.1 if r > 0.4 else 0.9
            dmax = min(max(dmax, 0.005), 0.3)
    log.append({"T": T, "E_per_N": n / len(s), "dmax": dmax, "t": time.time() - t0})
    print(json.dumps(log[-1]), flush=True)
    np.savez(f"{out}_T{T:.3f}.npz", s=s, box=box, lam=lam, T=T)
