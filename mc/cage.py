"""Local free-volume (cage) entropy of particles by environment.
Short NVT run from a snapshot (fixed box); per-particle positional covariance C_i over the run (unwrapped, COM removed).
s_i = 0.5*ln det(2*pi*e*C_i)  (Gaussian estimate of the cage entropy, k_B units).
Particles are classed by the tiles around them in the starting snapshot: 'tile' (only triangles/quads),
'defect' (touches a polygon with >=5 sides).  Also reports bond-length distributions by polygon type.
Usage: python cage.py snap.npz T nsweeps every"""
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "engine"))
import hcss_mc as mc
from analyze import tiles

fn, T = sys.argv[1], float(sys.argv[2]); nsw = int(sys.argv[3]); every = int(sys.argv[4])
z = np.load(fn); s = z["s"].copy(); box = np.asarray(z["box"], float); lam = float(z["lam"]); N = len(s)
X, ID, _, polys, _ = tiles(fn)
cls = np.zeros(N, int)  # 0 tile, 1 defect
bl = {3: [], 4: [], 5: [], 6: []}
for p in polys:
    n = len(p)
    if n >= 5:
        for v in p: cls[ID[v]] = 1
    key = n if n <= 6 else 6
    for i in range(n):
        bl[key].append(np.linalg.norm(X[p[(i + 1) % n]] - X[p[i]]))
img = np.zeros((N, 2), np.int64); n = mc.total_count(s, *box, lam); rs = np.random.default_rng(3)
traj = []; dmax = 0.02; acc = 0
for k in range(1, nsw + 1):
    n, ap, ab = mc.sweep_img(s, img, box, lam, 1 / T, 0.0, dmax, 0.0, 0.0, n, int(rs.integers(1, 2**31 - 1)))
    acc += ap
    if k % 200 == 0:
        dmax *= 1.1 if acc / (200 * N) > 0.4 else 0.9; acc = 0
    if k > nsw // 5 and k % every == 0:
        u = s + img; c = np.stack([box[0] * u[:, 0] + box[1] * u[:, 1], box[2] * u[:, 1]], 1)
        traj.append(c - c.mean(0))
traj = np.array(traj)  # (t, N, 2)
d = traj - traj.mean(0)
C = np.einsum("tni,tnj->nij", d, d) / len(traj)
det = np.linalg.det(C); msd = C[:, 0, 0] + C[:, 1, 1]
s_i = 0.5 * np.log((2 * math.pi * math.e) ** 2 * det)
out = dict(file=fn, N=N, frac_defect_particles=float(cls.mean()))
for name, m in (("tile", cls == 0), ("defect", cls == 1)):
    if m.any():
        out[name] = dict(n=int(m.sum()), msd=float(np.median(msd[m])), s_cage=float(np.mean(s_i[m])), s_cage_med=float(np.median(s_i[m])))
out["bond_len"] = {str(k): (float(np.mean(v)), float(np.percentile(v, 90))) for k, v in bl.items() if v}
out["max_msd"] = float(msd.max())
print(json.dumps(out, indent=1))
