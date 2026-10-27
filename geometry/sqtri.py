"""Square-triangle tilings: Stampfli inflation, faces, Q graph, diagonal cover."""
import numpy as np
from scipy.spatial import cKDTree
from collections import defaultdict, deque

S3 = np.sqrt(3.0)
LAM = 2 + S3
TOL = 1e-6


def dodecagon_filling(rot):
    """Vertices of a unit-edge dodecagon filled with 12 triangles + 6 squares.
    rot=0 or 1 (rotation by 30 deg)."""
    a0 = rot * np.pi / 6
    pts = [np.zeros(2)]
    hexv = [np.array([np.cos(a0 + k * np.pi / 3), np.sin(a0 + k * np.pi / 3)]) for k in range(6)]
    pts += hexv
    for k in range(6):
        p, q = hexv[k], hexv[(k + 1) % 6]
        n = (p + q) / 2
        n = n / np.linalg.norm(n)
        pts += [p + n, q + n]
    return np.array(pts)


DOD = [dodecagon_filling(0), dodecagon_filling(1)]


def dedupe(P, tol=1e-5):
    tree = cKDTree(P)
    pairs = tree.query_pairs(tol)
    parent = list(range(len(P)))

    def f(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j in pairs:
        a, b = f(i), f(j)
        if a != b:
            parent[max(a, b)] = min(a, b)
    keep = sorted({f(i) for i in range(len(P))})
    return P[keep]


def inflate(P, choose):
    """choose(i, p) -> 0/1 filling orientation for vertex i at position p."""
    out = [LAM * p + DOD[choose(i, p)] for i, p in enumerate(P)]
    return dedupe(np.vstack(out))


def build(P):
    """Unit-distance graph, rotation system, faces (triangles/squares)."""
    tree = cKDTree(P)
    pairs = np.array(sorted(tree.query_pairs(1 + 1e-4)))
    d = np.linalg.norm(P[pairs[:, 0]] - P[pairs[:, 1]], axis=1)
    pairs = pairs[np.abs(d - 1) < 1e-4]
    nbr = defaultdict(list)
    for i, j in pairs:
        nbr[i].append(j)
        nbr[j].append(i)
    # sort neighbours by angle (ccw)
    for i in nbr:
        v = P[nbr[i]] - P[i]
        ang = np.arctan2(v[:, 1], v[:, 0])
        nbr[i] = [nbr[i][k] for k in np.argsort(ang)]
    # faces: traverse half-edges, face on the left
    pos = {i: {j: k for k, j in enumerate(nbr[i])} for i in nbr}
    seen = set()
    faces = []
    for i in nbr:
        for j in nbr[i]:
            if (i, j) in seen:
                continue
            f = []
            a, b = i, j
            while (a, b) not in seen:
                seen.add((a, b))
                f.append(a)
                # next: at b, take the neighbour just clockwise from a  (face on left)
                k = pos[b][a]
                c = nbr[b][(k - 1) % len(nbr[b])]
                a, b = b, c
                if len(f) > 10:
                    break
            faces.append(f)
    def _area(f):
        x = P[f]
        return 0.5 * np.sum(x[:, 0] * np.roll(x[:, 1], -1) - np.roll(x[:, 0], -1) * x[:, 1])
    # drop clockwise 3-cycles: pieces of the outer face at boundary notches
    tris = [f for f in faces if len(f) == 3 and _area(f) > 0]
    sqs = [f for f in faces if len(f) == 4 and abs(np.linalg.norm(P[f[0]] - P[f[2]]) - np.sqrt(2)) < 1e-4]
    return nbr, tris, sqs


def interior_vertices(P, nbr, tris, sqs):
    ang = defaultdict(float)
    for f in tris:
        for v in f:
            ang[v] += np.pi / 3
    for f in sqs:
        for v in f:
            ang[v] += np.pi / 2
    return {v for v in nbr if abs(ang[v] - 2 * np.pi) < 1e-6}


def cover(P, nbr, sqs, interior):
    """Q graph and BFS construction of sigma. Returns sigma, conflicts, comps.
    Square s with corners f[0..3] (cyclic): c_s = +1 on f[0],f[2], -1 on f[1],f[3].
    sigma_s=+1 means diagonal f[0]-f[2]."""
    vsq = defaultdict(list)  # vertex -> [(s, c)]
    for s, f in enumerate(sqs):
        for k, v in enumerate(f):
            vsq[v].append((s, 1 if k % 2 == 0 else -1))
    edges = []  # (s, t, J, v)
    for v, L in vsq.items():
        if len(L) == 2:
            (s, cs), (t, ct) = L
            edges.append((s, t, -cs * ct, v))
    adj = defaultdict(list)
    for s, t, J, v in edges:
        adj[s].append((t, J, v))
        adj[t].append((s, J, v))
    sigma = {}
    comps = 0
    conflicts = []
    for s0 in range(len(sqs)):
        if s0 in sigma:
            continue
        comps += 1
        sigma[s0] = 1
        dq = deque([s0])
        while dq:
            s = dq.popleft()
            for t, J, v in adj[s]:
                if t not in sigma:
                    sigma[t] = sigma[s] * J
                    dq.append(t)
                elif sigma[t] != sigma[s] * J:
                    conflicts.append(v)
    return sigma, conflicts, comps, edges, vsq


def degrees_after(P, nbr, sqs, sigma):
    deg = {v: len(nbr[v]) for v in nbr}
    for s, f in enumerate(sqs):
        a, b = (f[0], f[2]) if sigma[s] == 1 else (f[1], f[3])
        deg[a] += 1
        deg[b] += 1
    return deg
