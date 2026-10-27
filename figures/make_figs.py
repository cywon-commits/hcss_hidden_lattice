import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib import font_manager as fm
from scipy.spatial import Voronoi, Delaunay, cKDTree
from collections import Counter, defaultdict
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "geometry"))
from sqtri import *

for f in fm.findSystemFonts():
    if "NotoSansCJK-Regular" in f:
        fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 15,
                     "axes.spines.top": False, "axes.spines.right": False})
OUT = HERE + "/"
import os
os.makedirs(OUT, exist_ok=True)

INK = "#1f2430"; MUTED = "#9aa0aa"; TRI = "#e9edf3"; SQ = "#c9dcf5"; ODD = "#e4572e"; EVEN = "#2d3142"
ACC = "#2a7f62"; WARM = "#f2a541"; BLUE = "#3a6ea5"


def save(fig, name):
    fig.savefig(OUT + name, dpi=200, bbox_inches="tight", pad_inches=0.05, transparent=True)
    plt.close(fig)


# ---------- f1: noisy triangular lattice Voronoi, odd cells ----------
def f1():
    rng = np.random.default_rng(7)
    a1, a2 = np.array([1, 0]), np.array([0.5, np.sqrt(3) / 2])
    pts = np.array([i * a1 + j * a2 for i in range(-14, 15) for j in range(-14, 15)])
    pts = pts[(np.abs(pts[:, 0]) < 9) & (np.abs(pts[:, 1]) < 7)]
    pts = pts + rng.normal(0, 0.21, pts.shape)
    vor = Voronoi(pts); tri = Delaunay(pts)
    deg = Counter()
    nb = defaultdict(set)
    for s in tri.simplices:
        for i in range(3):
            for j in range(3):
                if i != j: nb[s[i]].add(s[j])
    fig, ax = plt.subplots(figsize=(8, 6))
    win = lambda p: abs(p[0]) < 4.4 and abs(p[1]) < 2.9
    odd = []
    for i, r in enumerate(vor.point_region):
        reg = vor.regions[r]
        if -1 in reg or not reg or not win(pts[i]): continue
        z = len(nb[i])
        col = TRI if z == 6 else ("#f6c1ae" if z == 5 else "#b9cfe8") if z % 2 == 0 or True else TRI
        if z == 6: col = "#eef1f5"
        elif z == 5: col = "#f4b39c"
        elif z == 7: col = "#a9c4e6"
        else: col = "#d8d8d8"
        ax.add_patch(Polygon(vor.vertices[reg], closed=True, fc=col, ec="#6b7280", lw=0.8))
        if z % 2: odd.append(i)
    import networkx as nx
    G = nx.Graph()
    for i in odd:
        for j in odd:
            if i < j: G.add_edge(i, j, weight=-np.linalg.norm(pts[i] - pts[j]))
    for i, j in nx.max_weight_matching(G, maxcardinality=True):
        ax.plot(*zip(pts[i], pts[j]), color=ODD, lw=3.5, solid_capstyle="round", zorder=5)
    ax.set_xlim(-4.6, 4.6); ax.set_ylim(-3.1, 3.1); ax.set_aspect("equal"); ax.axis("off")
    save(fig, "f1_voronoi.png")


# ---------- square-triangle patch helpers ----------
rng = np.random.default_rng(3)
P2 = DOD[0].copy()
for k in range(2):
    P2 = inflate(P2, lambda i, p: int(rng.integers(2)))
nbr2, tris2, sqs2 = build(P2)
I2 = interior_vertices(P2, nbr2, tris2, sqs2)
sig2, conf2, comps2, edges2, vsq2 = cover(P2, nbr2, sqs2, I2)


def draw_tiling(ax, P, tris, sqs, box, sqcol=SQ, tricol=TRI, lw=0.7):
    for f in tris:
        if box(P[f].mean(0)): ax.add_patch(Polygon(P[f], fc=tricol, ec="#7a8190", lw=lw))
    for f in sqs:
        if box(P[f].mean(0)): ax.add_patch(Polygon(P[f], fc=sqcol, ec="#7a8190", lw=lw))


def f2():
    R = 5.2; box = lambda c: np.linalg.norm(c) < R
    fig, ax = plt.subplots(figsize=(7, 7))
    draw_tiling(ax, P2, tris2, sqs2, box)
    for v in nbr2:
        if np.linalg.norm(P2[v]) < R - 0.3 and v in I2:
            ax.plot(*P2[v], "o", ms=8, color=ODD if len(nbr2[v]) == 5 else EVEN, mec="white", mew=1)
    ax.set_xlim(-R + .4, R - .4); ax.set_ylim(-R + .4, R - .4); ax.set_aspect("equal"); ax.axis("off")
    save(fig, "f2_tiling.png")


def f3():
    R = 4.2; box = lambda c: np.linalg.norm(c) < R
    fig, ax = plt.subplots(figsize=(7, 7))
    draw_tiling(ax, P2, tris2, sqs2, box, sqcol="#dfe9f7")
    cen = np.array([P2[f].mean(0) for f in sqs2])
    for s, t, J, v in edges2:
        if box(cen[s]) and box(cen[t]):
            ax.plot([cen[s][0], P2[v][0], cen[t][0]], [cen[s][1], P2[v][1], cen[t][1]], color=MUTED, lw=1.2, ls=(0, (2, 2)))
    for s, f in enumerate(sqs2):
        if box(cen[s]):
            a, b = (f[0], f[2]) if sig2[s] == 1 else (f[1], f[3])
            ax.plot(*zip(P2[a], P2[b]), color=ACC, lw=3.2, solid_capstyle="round")
            ax.plot(*cen[s], "s", ms=6, color=BLUE)
    ax.set_xlim(-R + .5, R - .5); ax.set_ylim(-R + .5, R - .5); ax.set_aspect("equal"); ax.axis("off")
    save(fig, "f3_cover.png")


def develop_pos(P, tris, sqs, sigma, th):
    from develop import develop
    return develop(P, tris, sqs, sigma, th, return_positions=True)[0]


def f4():
    R = 3.6
    sel_t = [f for f in tris2 if np.linalg.norm(P2[f].mean(0)) < R]
    sel_s = [(s, f) for s, f in enumerate(sqs2) if np.linalg.norm(P2[f].mean(0)) < R]
    sq_l = [f for s, f in sel_s]; sg = {k: sig2[s] for k, (s, f) in enumerate(sel_s)}
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.2))
    for ax, th in zip(axs, (0, 30, 60)):
        Q = develop_pos(P2, sel_t, sq_l, sg, th)
        X = {v: Q[v] for v in Q}
        c0 = np.mean(list(X.values()), 0)
        for f in sel_t: ax.add_patch(Polygon([X[v] - c0 for v in f], fc=TRI, ec="#7a8190", lw=0.7))
        for k, f in enumerate(sq_l):
            col = SQ if th == 0 else ("#cfe8dc" if th == 30 else "#f7dcb0")
            ax.add_patch(Polygon([X[v] - c0 for v in f], fc=col, ec="#7a8190", lw=0.7))
            a, b = (f[0], f[2]) if sg[k] == 1 else (f[1], f[3])
            ax.plot(*zip(X[a] - c0, X[b] - c0), color=ACC, lw=1.6)
        allp = np.array([X[v] - c0 for v in X]); m = np.abs(allp).max() * 0.82
        ax.set_xlim(-m, m); ax.set_ylim(-m, m); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"θ = {th}°", fontsize=20, color=INK)
    save(fig, "f4_mechanism.png")


def f5():
    d = [0.5, 1.5, 3, 5, 7, 9, 11.5, 14.5, 18]
    data = {"k_d = 0.01": [1, .509, .266, .121, .06, .033, .017, .007, .003],
            "k_d = 0.1": [1, .287, .073, .009, 0, 0, .001, 0, 0],
            "k_d = 1": [1, .093, .007, 0, -.001, 0, 0, 0, 0]}
    cols = [ACC, BLUE, WARM]
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for (k, v), c in zip(data.items(), cols):
        ax.plot(d, v, "-o", color=c, lw=2.5, ms=6, label=k)
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xlabel("정사각형 사이 거리 (타일 변 단위)"); ax.set_ylabel("게이지 상관 ρ·g_s g_t")
    ax.legend(frameon=False); ax.set_xlim(0, 19)
    save(fig, "f5_corr.png")


def f6():
    fig, ax = plt.subplots(figsize=(7.5, 5))
    labels = ["이상 타일링", "무작위 대각선", "메커니즘에\n묶인 대각선"]
    vals = [4 * np.sqrt(3) - 6, 2 * np.sqrt(3) - 3, 0.0]
    b = ax.bar(labels, vals, color=[ODD, WARM, ACC], width=0.6)
    for r, v, t in zip(b, vals, ["0.928", "0.464", "≈ 0"]):
        ax.text(r.get_x() + r.get_width() / 2, v + 0.02, t, ha="center", fontsize=18, color=INK)
    ax.set_ylabel("홀수 배위 비율 n_odd / N"); ax.set_ylim(0, 1.05)
    save(fig, "f6_nodd.png")


def f7():
    LS = 2 * np.cos(np.radians(15))
    d = np.load(os.path.join(HERE, "..", "data", "configs", "lstar_approx2_sym.npz")); s = d["s"]; a, b, c = d["box"]
    M = np.array([[a, b], [0, c]]); P = s @ M.T; N = len(P)
    X = []; ID = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            X.append(P + i * M[:, 0] + j * M[:, 1]); ID.append(np.arange(N))
    X = np.vstack(X); ID = np.concatenate(ID)
    T = Delaunay(X).simplices
    corner = defaultdict(Counter)
    for sim in T:
        v = X[sim]; L = [np.linalg.norm(v[(k + 1) % 3] - v[(k + 2) % 3]) for k in range(3)]
        eq = all(abs(l - LS) < 1e-3 for l in L)
        for k in range(3):
            corner[ID[sim[k]]]["t" if eq else ("a" if abs(L[k] - 1) < 1e-3 else "o")] += 1
    ok = {i: corner[i]["a"] == corner[i]["o"] // 2 for i in range(N)}
    cen = M @ np.array([0.5, 0.5]); W = 9.5
    box = lambda p: abs(p[0] - cen[0]) < W and abs(p[1] - cen[1]) < W * 0.62
    fig, ax = plt.subplots(figsize=(9, 5.8))
    # rhombi: merge isosceles pairs across short (core) edge
    for sim in T:
        v = X[sim]
        if not box(v.mean(0)): continue
        L = [np.linalg.norm(v[(k + 1) % 3] - v[(k + 2) % 3]) for k in range(3)]
        eq = all(abs(l - LS) < 1e-3 for l in L)
        ax.add_patch(Polygon(v, fc=TRI if eq else "#f7dcb0", ec="#7a8190" if eq else "none", lw=0.6))
        if not eq:
            for k in range(3):
                if abs(L[k] - LS) < 1e-3:
                    p, q = v[(k + 1) % 3], v[(k + 2) % 3]; ax.plot(*zip(p, q), color="#7a8190", lw=0.6)
    for i in range(len(X)):
        if box(X[i]):
            ax.plot(*X[i], "o", ms=7, color=ACC if ok[ID[i]] else ODD, mec="white", mew=0.8)
    ax.set_xlim(cen[0] - W + .8, cen[0] + W - .8); ax.set_ylim(cen[1] - W * 0.62 + .8, cen[1] + W * 0.62 - .8)
    ax.set_aspect("equal"); ax.axis("off")
    save(fig, "f7_lstar.png")


def f8():
    """single 5-vertex: two sheared squares + 3 triangles, angles 90+θ / 90-θ"""
    th = np.radians(15)
    fig, ax = plt.subplots(figsize=(6, 6))
    v = np.zeros(2)
    seq = [("t", 60), ("t", 60), ("s+", 90 + 15), ("t", 60), ("s-", 90 - 15)]
    ang = 0
    for kind, deg in seq:
        a0 = np.radians(ang); a1 = np.radians(ang + deg)
        p = np.array([np.cos(a0), np.sin(a0)]); q = np.array([np.cos(a1), np.sin(a1)])
        if kind == "t":
            ax.add_patch(Polygon([v, p, q], fc=TRI, ec="#7a8190", lw=1.2))
        else:
            ax.add_patch(Polygon([v, p, p + q, q], fc=SQ, ec="#7a8190", lw=1.2))
            if kind == "s+": ax.plot(*zip(v, p + q), color=ACC, lw=3.5)
            else: ax.plot(*zip(p, q), color=ACC, lw=3.5)
            mid = np.radians(ang + deg / 2)
            rr, dd = (0.6, 0.3) if kind == 's+' else (0.42, -0.05)
            ax.text(rr * np.cos(mid+dd), rr * np.sin(mid+dd), f"{90 + 15 if kind == 's+' else 90 - 15}°" if False else ("90°+θ" if kind == "s+" else "90°−θ"),
                    ha="center", va="center", fontsize=17, color=INK, fontweight="bold")
        if kind == "t":
            mid = np.radians(ang + deg / 2)
            ax.text(0.45 * np.cos(mid), 0.45 * np.sin(mid), "60°", ha="center", va="center", fontsize=14, color=MUTED)
        ang += deg
    ax.plot(0, 0, "o", ms=12, color=ODD, mec="white", mew=1.5)
    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6); ax.set_aspect("equal"); ax.axis("off")
    save(fig, "f8_vertex.png")


if __name__ == "__main__":
    for fn in (f1, f2, f3, f4, f5, f6, f7, f8):
        fn(); print("done", fn.__name__)
