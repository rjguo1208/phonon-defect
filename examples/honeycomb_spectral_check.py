#!/usr/bin/env python3
"""Brute-force check of the honeycomb spectral function at a finite concentration (Honeycomb page).

Random substitutions on A sites, n_d = 0.02 per cell, on a 36 x 36 torus (26 defects): every
configuration is diagonalized exactly, its displacement eigenvectors are unfolded onto Bloch waves of
the host (host-mass weighted), and A(q, w) is averaged over 10 configurations and over the star of q.
Compared with the dilute-limit forms of scripts/honeycomb_figures.py (eta = 0.02): the Dyson equation,
the first-order average, and the rule used for the figure (Dyson where the host has modes, first order
where it has none). About 12 minutes on 16 cores.

Run:  python examples/honeycomb_spectral_check.py
"""
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from phdef import honeycomb as hc  # noqa: E402
import honeycomb_figures as hf  # noqa: E402

L, NCONF, ETA, ND, LQ = 36, 10, 0.02, 0.02, 256
OM = np.arange(0.005, 3.4, 0.005)
BMAT = np.c_[hc.B1, hc.B2]


def path_points():
    """Grid points (m1, m2), q = (m1 b1 + m2 b2) / L, along Γ-M-K-Γ (L divisible by 6)."""
    g, m, k = np.array([0, 0]), np.array([L // 2, 0]), np.array([2 * L // 3, L // 3])
    seg = lambda a, b, n: [tuple(np.round(a + (b - a) * t / n).astype(int)) for t in range(n)]
    return seg(g, m, L // 2) + seg(m, k, L // 6) + seg(k, g, L // 3) + [(0, 0)]


def star(m):
    q = BMAT @ np.array(m, float) / L
    out = set()
    for n in range(6):
        c, s = np.cos(np.pi * n / 3), np.sin(np.pi * n / 3)
        for mirror in (1, -1):
            r = np.linalg.solve(BMAT, np.array([[c, -s], [s, c]]) @ np.array([q[0], mirror * q[1]])) * L
            if np.abs(r - np.round(r)).max() < 1e-6:
                out.add((int(round(r[0])) % L, int(round(r[1])) % L))
    return sorted(out)


def configuration(Md, fac, seed, stars):
    rng = np.random.default_rng(seed)
    isdef = np.zeros(2 * L * L, bool)
    isdef[2 * rng.choice(L * L, round(ND * L * L), replace=False)] = True     # A atoms
    n1, n2 = (a.ravel() for a in np.meshgrid(np.arange(L), np.arange(L), indexing="ij"))
    rows, cols, vals = [], [], []
    for atoms, C in hc.terms():
        ids = [hc.site_index(s, n1 + c[0], n2 + c[1], L) for s, c in atoms]
        factor = np.where(np.any(isdef[np.array(ids)], axis=0), fac, 1.0)      # every term with a defect atom
        for a in range(len(atoms)):
            for b in range(len(atoms)):
                for al, be in zip(*np.nonzero(C[a][b])):
                    rows.append(3 * ids[a] + al)
                    cols.append(3 * ids[b] + be)
                    vals.append(factor * C[a][b][al, be])
    Phi = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(6 * L * L,) * 2).toarray()
    m = np.tile([1.0, 1.3], L * L)
    m[isdef] = Md
    w = np.repeat(1.0 / np.sqrt(m), 3)
    lam, X = np.linalg.eigh(Phi * np.outer(w, w))
    Y = (X * (w * np.repeat(np.sqrt(np.tile([1.0, 1.3], L * L)), 3))[:, None]).reshape(L, L, 6, -1)
    W = (np.abs(np.fft.fft2(Y, axes=(0, 1)) / L) ** 2).sum(axis=2)          # weight of each mode at each q
    Wp = np.array([np.mean([W[a, b] for a, b in st], axis=0) for st in stars])
    return Wp @ (-(2 * OM[None, :] / np.pi) * (1.0 / ((OM[None, :] + 1j * ETA) ** 2 - lam[:, None])).imag)


def main():
    t0 = time.time()
    path = path_points()
    stars = [star(m) for m in path]
    qs = np.array([BMAT @ np.array(m, float) / L for m in path])
    labels = {"Γ": 0, "M": L // 2, "K": L // 2 + L // 6}
    offs, subl, _, _ = hc.perturbation((0, (0, 0), 0.3, 1.3))
    g0 = hc.host_green_support(OM, ETA, LQ, offs, subl)
    near = (OM > 1.30) & (OM < 1.60)
    acoustic = OM < 0.5                    # the flexural branch, w^2 ~ q^4, is the most sensitive to the phase of Sigma(q)
    print(f"{L} x {L} torus, {round(ND * L * L)} A-site defects, {NCONF} configurations, eta = {ETA}")
    for key, _, Md, fac in hf.DEFECTS:
        brute = np.mean([configuration(Md, fac, 100 * (key == "heavy") + c, stars) for c in range(NCONF)], axis=0)
        offs, subl, dPhi, dM = hc.perturbation((0, (0, 0), Md, fac))
        T = hc.tmatrix_u(OM, ETA, g0, dPhi, dM)
        rule = sum(hf.averaged_spectral(qs, OM, ETA, T, offs, subl, ND))
        saved = hf.FREE
        hf.FREE = []                                                        # Dyson everywhere
        dyson = sum(hf.averaged_spectral(qs, OM, ETA, T, offs, subl, ND))
        hf.FREE = [(0.0, np.inf)]                                           # first order everywhere
        first = sum(hf.averaged_spectral(qs, OM, ETA, T, offs, subl, ND))
        hf.FREE = saved
        print(f"\n== {key} ({time.time() - t0:.0f} s)")
        for name, A in (("Dyson", dyson), ("first order", first), ("rule", rule)):
            l1 = lambda w: np.abs(brute[:, w] - A[:, w]).sum(1) / np.abs(brute[:, w]).sum(1)
            ridge = (dyson > 1.0) & ((OM[None, :] < 1.36) | ((OM[None, :] > 1.55) & (OM[None, :] < 1.86)))
            print(f"   {name:11s}: L1 distance in 1.30-1.60, median over the path {np.median(l1(near)):.3f}; "
                  f"below 0.5 (acoustic) median {np.median(l1(acoustic)):.3f}, max {l1(acoustic).max():.3f}; "
                  f"host ridge median deviation {np.median((np.abs(brute - A) / dyson)[ridge]):.3f}")
        for lo, hi in ((1.42, 1.50), (2.6, 3.4)):
            sel = (OM > lo) & (OM < hi)
            for lab, i in labels.items():
                pk = lambda A: f"{OM[sel][np.argmax(A[i][sel])]:.3f}/{np.trapezoid(A[i][sel], OM[sel]):.4f}"
                print(f"   peak/weight in {lo}-{hi} at {lab}: brute {pk(brute)}, Dyson {pk(dyson)}, rule {pk(rule)}")


if __name__ == "__main__":
    main()
