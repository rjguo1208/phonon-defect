#!/usr/bin/env python3
"""Brute-force check of the spectral function at a finite defect concentration (Toy page, Figure 5).

Random defects at n_d = 0.02 on a 64 x 64 torus: every configuration is diagonalized exactly, its
displacement eigenvectors are unfolded onto plane waves, and A(k, w) is averaged over 12 configurations
and over the star of each k. The result is compared with the dilute-limit A(k, w) of
scripts/toy_figures.py (Dyson equation with Sigma = n_d T in the host band, first-order average above
the middle of the gap), and with the Dyson equation used everywhere.

Run:  python examples/toy_spectral_check.py      (about two minutes on 8 cores)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from toy_figures import host_green, defect_tmatrix, spectral_function  # noqa: E402

L, NCONF, ETA = 64, 12, 0.01
ND = round(0.02 * L * L) / (L * L)                         # 82 defects
DEFECTS = {"light": (0.3, 1.3), "heavy": (8.0, 0.6)}
W_CUT = 0.5 * (2 * np.sqrt(2) + 4.348826)                  # middle of the gap below the local mode
OM = np.arange(0.002, 4.9, 0.004)


def path_indices():
    """Grid points (p, q), k = 2 pi (p, q) / L, along Γ-X-M-Γ."""
    h = L // 2
    return [(n, 0) for n in range(h)] + [(h, n) for n in range(h)] + [(n, n) for n in range(h, -1, -1)]


def star(p, q):
    return sorted({((sa * a) % L, (sb * b) % L) for a, b in ((p, q), (q, p)) for sa in (1, -1) for sb in (1, -1)})


def configuration(Md, Kd, seed, path):
    """A(k, w) of one random configuration, averaged over the star of each path point."""
    N = L * L
    rng = np.random.default_rng(seed)
    isdef = np.zeros(N, bool)
    isdef[rng.choice(N, round(ND * N), replace=False)] = True
    i, j = np.divmod(np.arange(N), L)
    Phi = np.zeros((N, N))
    a = np.arange(N)
    for di, dj in ((1, 0), (0, 1)):
        b = ((i + di) % L) * L + (j + dj) % L
        k = np.where(isdef | isdef[b], Kd, 1.0)            # a bond touching a defect has spring K'
        Phi[a, b] -= k
        Phi[b, a] -= k
        Phi[a, a] += k
        Phi[b, b] += k
    s = 1.0 / np.sqrt(np.where(isdef, Md, 1.0))
    lam, X = np.linalg.eigh(s[:, None] * Phi * s[None, :])
    W = np.abs(np.fft.fft2((s[:, None] * X).reshape(L, L, N), axes=(0, 1)) / L) ** 2   # |<k|u_beta>|^2
    Wp = np.array([np.mean([W[pq] for pq in star(*pq0)], axis=0) for pq0 in path])
    kernel = -(2 * OM[None, :] / np.pi) * (1.0 / ((OM[None, :] + 1j * ETA) ** 2 - lam[:, None])).imag
    return Wp @ kernel


def peaks(a, frac):
    return [(round(float(OM[j]), 3), round(float(a[j]), 1)) for j in range(1, len(a) - 1)
            if a[j] > a[j - 1] and a[j] > a[j + 1] and a[j] > frac * a.max()]


def main():
    t0 = time.time()
    path = path_indices()
    kvec = 2 * np.pi * np.array(path, float) / L
    h = L // 2
    g = host_green(OM, ETA)
    print(f"{L} x {L} torus, {round(ND * L * L)} defects (n_d = {ND:.5f}), {NCONF} configurations, eta = {ETA}")
    for name, (Md, Kd) in DEFECTS.items():
        brute = np.mean([configuration(Md, Kd, 1000 * (name == "heavy") + c, path) for c in range(NCONF)], axis=0)
        T, _ = defect_tmatrix(OM, ETA, Md, Kd, g)
        dilute = spectral_function(kvec, OM, ETA, T, ND, W_CUT)
        dyson = spectral_function(kvec, OM, ETA, T, ND)
        band = (dilute > 1.0) & (OM[None, :] < 2.9)
        print(f"\n== {name} defect ({time.time() - t0:.0f} s): median relative deviation on the host ridge "
              f"{np.median((np.abs(brute - dilute) / dilute)[band]):.3f}")
        if name == "light":
            sel = OM > 3.5
            print("   local-mode band, peak and weight above w = 3.5:  brute force | dilute | Dyson everywhere")
            for lab, i in (("Γ", 0), ("X", h), ("M", 2 * h)):
                row = [(OM[sel][np.argmax(A[i][sel])], np.trapezoid(A[i][sel], OM[sel])) for A in (brute, dilute, dyson)]
                print(f"   {lab}: " + " | ".join(f"{p:.3f} / {w:.4f}" for p, w in row))
        else:
            print("   peaks (w, A) along Γ-X below w = 1.2:  brute force | dilute")
            sel = OM < 1.2
            for n in range(2, 11):
                print(f"   k = {kvec[n, 0]:.3f}: {peaks(brute[n][sel], 0.15)} | {peaks(dilute[n][sel], 0.15)}")


if __name__ == "__main__":
    main()
