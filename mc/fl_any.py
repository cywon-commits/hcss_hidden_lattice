"""Frenkel–Ladd (Einstein molecule, Vega & Noya 2007) Gibbs free energy per particle of the basin of a given
equilibrated HCSS configuration (npz with s, box, lam).  Same formula as hcss_lambda_star/sim/fl_free_energy.py,
but sites = NPT-averaged (unwrapped, drift-removed) positions of the given configuration, no rescaling.
Usage: python fl_any.py config.npz T P out.json [npt_sweeps] [sample]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "engine"))
import hcss_mc as mc


def block_err(x, nb=10):
    x = np.asarray(x, float); m = len(x) // nb
    if m < 1: return float("nan")
    b = x[: m * nb].reshape(nb, m).mean(1); return float(b.std(ddof=1) / math.sqrt(nb))


def main():
    fn, T, P, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    nnpt = int(sys.argv[5]) if len(sys.argv) > 5 else 40000
    nsamp = int(sys.argv[6]) if len(sys.argv) > 6 else 6000
    z = np.load(fn); s = z["s"].copy(); box = np.asarray(z["box"], float).copy(); lam = float(z["lam"])
    beta = 1 / T; N = len(s); rs = np.random.default_rng(7); t0 = time.time()
    img = np.zeros((N, 2), np.int64); n = mc.total_count(s, *box, lam); assert n >= 0
    dmax, dbox = 0.02, 1e-4; acc_p = acc_b = win = 0
    boxes, es, vs = [], [], []; u_acc = np.zeros_like(s); u_cnt = 0; snaps = []
    for k in range(nnpt):
        n, ap_, ab_ = mc.sweep_img(s, img, box, lam, beta, P, dmax, dbox, 0.0, n, int(rs.integers(1, 2**31 - 1)))
        acc_p += ap_; acc_b += ab_; win += 1
        if k < nnpt // 4 and win == 100:
            dmax *= 1.1 if acc_p / (win * N) > 0.35 else 0.9
            dbox *= 1.1 if acc_b / win > 0.3 else 0.9
            acc_p = acc_b = win = 0
        if k >= nnpt // 2:
            boxes.append(box.copy()); es.append(n / N); vs.append(box[0] * box[2] / N)
            if k % 10 == 0:
                u = s + img; u = u - u.mean(0); u_acc += u; u_cnt += 1
                if k % 200 == 0: snaps.append(u.copy())
    boxm = np.mean(boxes, 0); e_m = float(np.mean(es)); v_m = float(np.mean(vs))
    um = u_acc / u_cnt; s0 = (um + 0.5) % 1.0
    nc = mc.total_count(s0, *boxm, lam); used = "mean"
    if nc < 0:
        d = [((((x - um) + 0.5) % 1.0 - 0.5) ** 2).sum() for x in snaps]
        s0 = (snaps[int(np.argmin(d))] + 0.5) % 1.0; used = "medoid"
        if mc.total_count(s0, *boxm, lam) < 0:
            # last configuration rescaled into the averaged box (fractional coordinates are box-independent)
            s0 = s.copy(); used = "last"
            if mc.total_count(s0, *boxm, lam) < 0:
                boxm = box.copy(); V = None; used = "last+box"

    V = boxm[0] * boxm[2]; assert mc.total_count(s0, *boxm, lam) >= 0; lam_min = 0.05; lam_max = 3.0e4
    for _ in range(8):
        dA1, free, _m = mc.einstein_dA1(s0, boxm, lam, beta, lam_max, 2000, 3)
        if free >= 0.8: break
        lam_max *= 3.0
    xg, wg = np.polynomial.legendre.leggauss(14)
    x0, x1 = math.log(lam_min), math.log(lam_max)
    xs = 0.5 * (x1 - x0) * xg + 0.5 * (x1 + x0); ws = 0.5 * (x1 - x0) * wg
    lams = [float(v) for v in np.exp(xs)] + [lam_min]
    msd, msd_err = {}, {}
    for j in sorted(range(len(lams)), key=lambda j: -lams[j]):
        L = lams[j]; s = s0.copy(); n = mc.total_count(s, *boxm, lam)
        dm = min(0.05, 0.6 / math.sqrt(L)); acc = cnt = 0; vals = []
        for k in range(2000 + nsamp):
            n, ac, tot = mc.sweep_fl(s, s0, boxm, lam, beta, L, dm, n, int(rs.integers(1, 2**31 - 1)))
            if k < 2000:
                acc += ac; cnt += 1
                if cnt == 100:
                    dm *= 1.1 if acc / (cnt * (N - 1)) > 0.4 else 0.9; acc = cnt = 0
            else:
                vals.append(tot)
        msd[j] = float(np.mean(vals)); msd_err[j] = block_err(vals)
    K = len(xs)
    I = sum(ws[j] * lams[j] * msd[j] for j in range(K)) + lam_min * msd[K]
    I_err = math.sqrt(sum((ws[j] * lams[j] * msd_err[j]) ** 2 for j in range(K)) + (lam_min * msd_err[K]) ** 2)
    bA = -math.log(V) - (N - 1) * math.log(math.pi / lam_max) + dA1 - I
    bg = bA / N + beta * P * v_m
    res = dict(config=fn, T=T, P=P, N=N, e=e_m, v=v_m, rho=1 / v_m, h=e_m + P * v_m, sites=used,
               beta_g=bg, beta_g_err=math.hypot(I_err / N, beta * P * block_err(vs)), lam_max=lam_max,
               dA1_free=free, msd_low=msd[K], seconds=round(time.time() - t0))
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("config", "P", "rho", "e", "h", "beta_g", "beta_g_err", "sites", "msd_low", "seconds")}), flush=True)


if __name__ == "__main__":
    main()
