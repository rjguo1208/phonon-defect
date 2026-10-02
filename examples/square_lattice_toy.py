#!/usr/bin/env python3
"""Toy validation of the defect-phonon upfolding (cluster + Krylov chain of the host).

2D square lattice, scalar displacements, K = 1, M = 1, L x L = 121 x 121 sites.
Defect at the centre: mass M' and spring K' on its four bonds. Cluster C = Manhattan
ball of radius 2 (13 sites) -- contains every mass/spring change, so the coupling
block D_BC and the bulk D_BB are host quantities and the chain is built ONCE.

Checks, for a light defect (local mode above the band) and a heavy defect
(low-frequency resonance inside the band):
  * top mode of Ct vs the exact top eigenvalue of the full lattice;
  * G_00(omega + i eta) from the dressed-mode sum vs a direct sparse solve;
  * mode sum vs block continued fraction (same chain) -- identical by construction.

Run:  python examples/square_lattice_toy.py      (~1 min, mostly the reorthogonalised chain)
"""
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from phdef import block_lanczos, assemble, modes, cluster_green_modes, cluster_green_cf  # noqa: E402
from phdef.models import square_lattice_dynmat, manhattan_partition  # noqa: E402

L, RC, MMAX, ETA = 121, 2, 160, 0.05
OMEGAS = np.array([0.15, 0.25, 0.35, 0.5, 0.8, 1.2, 1.8, 2.5])
CASES = [("light defect M'=0.3, K'=1.3 (local mode above the band)", 0.3, 1.3),
         ("heavy defect M'=8,   K'=0.6 (resonance inside the band)", 8.0, 0.6)]


def main():
    C, Bk = manhattan_partition(L, RC)
    host, c0 = square_lattice_dynmat(L)
    i0 = int(np.where(C == c0)[0][0])
    D_BB = host[Bk][:, Bk].tocsr()
    D_BC = host[Bk][:, C].toarray()

    t0 = time.time()
    B0, A_all, B_all = block_lanczos(D_BB, D_BC, MMAX)
    print(f"host chain: |C|={len(C)}, |B|={len(Bk)}, block width {B0.shape[0]}, "
          f"{len(A_all)} blocks, built once in {time.time() - t0:.1f} s")

    z = (OMEGAS + 1j * ETA) ** 2
    for name, Md, Kd in CASES:
        D, _ = square_lattice_dynmat(L, Md, Kd)
        assert np.allclose(D[Bk][:, C].toarray(), D_BC), "coupling block must be defect independent"
        D_CC = D[C][:, C].toarray()
        top_exact = sla.eigsh(D, k=1, which="LA", return_eigenvectors=False)[0]
        e0 = np.zeros(D.shape[0], dtype=complex)
        e0[c0] = 1.0
        Dcsc = D.tocsc()
        I = sp.identity(D.shape[0], format="csc")
        G_exact = np.array([sla.spsolve(zz * I - Dcsc, e0)[c0] for zz in z])
        print(f"\n== {name}")
        for m in (10, 40, MMAX):
            Ct = assemble(D_CC, B0, A_all[:m], B_all[:m - 1])
            w2, V = modes(Ct)
            G_modes = cluster_green_modes(w2, V, len(C), z)[:, i0, i0]
            err = np.max(np.abs(G_modes - G_exact) / np.abs(G_exact))
            line = (f"   m={m:3d} (dim {Ct.shape[0]:5d}): top mode {np.sqrt(w2[-1]):.6f} "
                    f"(exact {np.sqrt(top_exact):.6f}), max rel err G_00 = {err:.1e}")
            if m == 40:
                G_cf = cluster_green_cf(D_CC, B0, A_all[:m], B_all[:m - 1], z)[:, i0, i0]
                line += f", |CF - modes| = {np.max(np.abs(G_cf - G_modes)):.1e}"
            print(line)
        rho = -2 * OMEGAS / np.pi * G_exact.imag
        print("   exact local DOS rho_00:", " ".join(f"{w:.2f}:{x:.3f}" for w, x in zip(OMEGAS, rho)))
    print(f"\nhost band top sqrt(8) = {np.sqrt(8):.6f}")


if __name__ == "__main__":
    main()
