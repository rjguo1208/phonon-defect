#!/usr/bin/env python3
"""Figures and data for the honeycomb toy page (run: npm run honeycomb).

--host     site/results/honeycomb-dispersion.* : host dispersion beside the density of states,
           split into in-plane and out-of-plane motion (seconds).
--defects  light and heavy substitution: chain of the bath from the Bloch modes of a 64 x 64
           q-mesh, upfolded defect Green functions against the exact 64 x 64 supercell; writes
           the numbers as CSV under site/data/ (a few minutes; run it on a compute node).
--plot     site/results/honeycomb-ldos.* and honeycomb-convergence.* from those CSV files.
Without options these three run.
--polar    the same chain on the polar host (32 x 32 mesh, minutes), Figures 5 and 6;
           --polar-plot redraws them from the CSV file.
--spectral site/results/honeycomb-spectral.* : configuration-averaged A(q, w) at n_d = 0.02 per
           cell from the T-matrix of one substitution (hc.perturbation, hc.tmatrix_u,
           hc.self_energy; 15 minutes on 8 cores); --spectral-plot redraws it from the .npz file
           named by HONEYCOMB_SPECTRAL_CACHE.
Same style and colours as scripts/toy_figures.py.
"""
import argparse
import csv
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import matplotlib

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from phdef import honeycomb as hc  # noqa: E402
from phdef import assemble, modes, cluster_green_cf  # noqa: E402
from toy_figures import plt, INK, INK2, MUTED, AXIS, save, write_csv  # noqa: E402

COLOR = {"in": "#2a78d6", "out": "#eb6834"}          # categorical slots 1-2, validated
NGRID, SIGMA = 300, 0.006


def branches(k):
    """In-plane (4) and out-of-plane (2) frequencies at the k-points k (n, 2)."""
    D = np.array([hc.bloch(kk) for kk in k])
    w = lambda idx: np.sqrt(np.clip(np.linalg.eigvalsh(D[:, idx][:, :, idx]), 0.0, None))
    return w(hc.IN_PLANE), w(hc.OUT_OF_PLANE)


def host_figure():
    s, kpath, ticks = hc.path(per_unit=60)
    w_in, w_out = branches(kpath)
    i, j = np.meshgrid(np.arange(NGRID), np.arange(NGRID), indexing="ij")
    kgrid = (np.outer(i.ravel(), hc.B1) + np.outer(j.ravel(), hc.B2)) / NGRID
    g_in, g_out = branches(kgrid)
    w_dos = np.linspace(0.0, 2.0, 1601)
    dw = w_dos[1] - w_dos[0]
    edges = np.r_[w_dos - dw / 2, w_dos[-1] + dw / 2]
    x = np.arange(-5 * SIGMA, 5 * SIGMA + dw / 2, dw)
    kernel = np.exp(-0.5 * (x / SIGMA) ** 2)
    kernel /= kernel.sum()
    dos = lambda w: np.convolve(np.histogram(w.ravel(), bins=edges)[0] / (len(kgrid) * dw), kernel, mode="same")
    rho_in, rho_out = dos(g_in), dos(g_out)
    top = max(g_in.max(), g_out.max())

    fig, (ax, axd) = plt.subplots(1, 2, figsize=(10.0, 4.4), sharey=True, constrained_layout=True,
                                  gridspec_kw=dict(width_ratios=[1.85, 1.0]))
    for n in range(4):
        ax.plot(s, w_in[:, n], color=COLOR["in"], lw=1.5, solid_capstyle="round", label="in-plane (x, y)" if n == 0 else None)
    for n in range(2):
        ax.plot(s, w_out[:, n], color=COLOR["out"], lw=1.5, solid_capstyle="round", label="out-of-plane (z)" if n == 0 else None)
    ax.set_xlim(0, ticks[-1])
    ax.set_xticks(ticks, ["Γ", "M", "K", "Γ"])
    ax.grid(axis="x", visible=False)
    for t in ticks[1:-1]:
        ax.axvline(t, color=AXIS, lw=0.75)
    ax.set_ylim(0, 2.0)
    ax.set_ylabel("Frequency ω")
    ax.set_title("Host phonon dispersion", loc="left", color=INK, fontsize=11)
    ax.legend(loc="lower center", ncol=2, fontsize=9.5, labelcolor=INK2, bbox_to_anchor=(0.48, 0.0))
    x0 = ticks[1] - 0.2                                           # direct labels just before M
    j0 = int(np.argmin(np.abs(s - x0)))
    for lab, y in (("ZA", w_out[j0, 0]), ("TA", w_in[j0, 0]), ("ZO", w_out[j0, 1]), ("LA", w_in[j0, 1])):
        ax.text(x0, y + 0.025, lab, color=INK2, fontsize=9, ha="right", va="bottom")
    ax.text(0.05, w_in[0, 3] + 0.03, "LO, TO", color=INK2, fontsize=9, ha="left", va="bottom")
    ins = ax.inset_axes([0.035, 0.58, 0.15, 0.18])                # Brillouin zone and the path
    ins.set_facecolor("none")
    corners = np.array([[np.cos(np.pi / 3 * n), np.sin(np.pi / 3 * n)] for n in range(7)]) * np.linalg.norm(hc.POINTS["K"])
    ins.plot(corners[:, 0], corners[:, 1], color=AXIS, lw=0.75)
    pk = np.array([hc.POINTS[p] for p in ("Γ", "M", "K", "Γ")])
    ins.plot(pk[:, 0], pk[:, 1], color=INK, lw=1.2, solid_joinstyle="round")
    for lab, (x, y), ha, va in (("Γ", pk[0], "right", "bottom"), ("M", pk[1], "left", "top"), ("K", pk[2], "left", "center")):
        ins.text(x + (0.15 if ha == "left" else -0.15), y, lab, ha=ha, va=va, color=INK2, fontsize=9)
    ins.set_aspect("equal")
    ins.set_xlim(-3.2, 3.6)
    ins.set_ylim(-3.0, 3.0)
    ins.axis("off")

    axd.plot(rho_in, w_dos, color=COLOR["in"], lw=1.5, label="in-plane")
    axd.plot(rho_out, w_dos, color=COLOR["out"], lw=1.5, label="out-of-plane")
    axd.set_xlim(0, None)
    axd.set_xlabel("Density of states per cell")
    axd.set_title("Density of states", loc="left", color=INK, fontsize=11)
    axd.tick_params(labelleft=False)
    axd.legend(loc="upper right", fontsize=9.5, labelcolor=INK2, bbox_to_anchor=(1.0, 0.63))
    gap = (float(g_in[:, 1].max()), float(g_in[:, 2].min()))       # LA top to the bottom of the optical branches
    axd.text(3.0, 0.5 * sum(gap), f"gap {gap[0]:.2f}–{gap[1]:.2f}", color=INK2, fontsize=9, ha="left", va="center")
    for a in (ax, axd):
        a.axhline(top, color=AXIS, lw=0.75)
    save(fig, "honeycomb-dispersion")

    write_csv("honeycomb-dispersion.csv", ["path_coordinate", "kx", "ky"] + [f"omega_in_{n + 1}" for n in range(4)]
              + [f"omega_out_{n + 1}" for n in range(2)],
              [(f"{a:.6f}", f"{kk[0]:.6f}", f"{kk[1]:.6f}", *[f"{v:.6f}" for v in (*wi, *wo)])
               for a, kk, wi, wo in zip(s, kpath, w_in, w_out)])
    write_csv("honeycomb-dos.csv", ["omega", "rho_in_plane", "rho_out_of_plane"],
              [(f"{w:.5f}", f"{a:.6e}", f"{b:.6e}") for w, a, b in zip(w_dos, rho_in, rho_out)])

    for lab in ("Γ", "M", "K"):
        wi, wo = branches(hc.POINTS[lab][None, :])
        print(f"{lab}: in-plane {np.round(wi[0], 4)}, out-of-plane {np.round(wo[0], 4)}")
    print("branch ranges on the grid:", [(round(float(g.min()), 4), round(float(g.max()), 4))
                                         for g in (*g_in.T, *g_out.T)])
    print(f"top of the spectrum {top:.4f}; DOS weights in/out = {np.trapezoid(rho_in, w_dos):.4f} / "
          f"{np.trapezoid(rho_out, w_dos):.4f}")


# ---------------------------------------------------------------- defects
L_MESH, MMAX, ETA_D, RADIUS = 64, 240, 0.02, 2
DEFECTS = [("light", "Light defect", 0.3, 1.3), ("heavy", "Heavy defect", 8.0, 0.6)]
TEST = np.array([0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9, 1.2, 1.45, 1.6, 1.75, 2.5])
M_SCAN = [5, 10, 20, 40, 60, 80, 120, 160, 200, 240]
M_LOCAL = list(range(1, 13)) + [20, 40, 80]
LOCAL = [("in-plane local mode", 3.155, 2), ("out-of-plane mode in the gap", 1.4543, 1),
         ("in-plane mode in the gap", 1.3773, 1)]            # light defect: targets near the exact values
DATA = os.path.join(ROOT, "site", "data")


def site_green(D, coords, omegas, eta):
    """Exact (z - D)^-1 between the given coordinates, one sparse LU per frequency."""
    n = D.shape[0]
    A0 = D.tocsc().astype(complex)
    I = sp.identity(n, format="csc", dtype=complex)
    rhs = np.zeros((n, len(coords)), complex)
    rhs[coords, np.arange(len(coords))] = 1.0
    return np.array([sla.splu((((w + 1j * eta) ** 2) * I - A0).tocsc()).solve(rhs)[coords] for w in omegas])


def rho(omega, G):
    """Local density of states per coordinate from the diagonal of G: (x, z)."""
    return (-(2 * omega / np.pi) * G[:, 0, 0].imag, -(2 * omega / np.pi) * G[:, 2, 2].imag)


def defect_data():
    t0 = time.time()
    L = L_MESH
    site = hc.site_index(0, 0, 0, L)
    hop = hc.hop_distance(site, L)
    C, Bk = hc.coordinates(np.where(hop <= RADIUS)[0]), hc.coordinates(np.where(hop > RADIUS)[0])
    own = [int(np.where(C == 3 * site + a)[0][0]) for a in range(3)]          # defect coordinates in C
    cache = os.environ.get("HONEYCOMB_CHAIN_CACHE")          # optional .npz outside the repository
    if cache and os.path.exists(cache):
        z = np.load(cache)
        B0, A1, W = z["B0"], z["A1"], z["W"]
        A, Bc = [z[f"A{j}"] for j in range(MMAX)], [z[f"B{j}"] for j in range(MMAX - 1)]
    else:
        B0, A, Bc, (A1, W) = hc.mode_chain(C, L, MMAX)
        if cache:
            np.savez(cache, B0=B0, A1=A1, W=W, **{f"A{j}": a for j, a in enumerate(A)}, **{f"B{j}": b for j, b in enumerate(Bc)})
    D0, _ = hc.supercell(L)
    print(f"chain from the Bloch modes of the {L} x {L} mesh: {len(A)} blocks, widths {A[0].shape[0]}..{A[-1].shape[0]}, "
          f"first block check {np.abs(A1 - W.conj().T @ D0[C][:, C].toarray() @ W).max():.1e}, {time.time() - t0:.0f} s")
    conv, local, ldos_up, ldos_ex = [], [], [], []
    w_up = np.linspace(0.005, 3.3, 1320)
    for key, _, Md, fac in DEFECTS:
        D, _ = hc.supercell(L, defect=(0, (0, 0), Md, fac))
        assert abs(D[Bk][:, C] - D0[Bk][:, C]).max() < 1e-14
        Dcc = D[C][:, C].toarray()
        Gx = site_green(D, 3 * site + np.arange(3), TEST, ETA_D)
        for m in M_SCAN:
            G = cluster_green_cf(Dcc, B0, A[:m], Bc[:m - 1], (TEST + 1j * ETA_D) ** 2)[:, own][:, :, own]
            err_in = max(np.max(np.abs(G[:, a, a] - Gx[:, a, a]) / np.abs(Gx[:, a, a])) for a in (0, 1))
            err_out = np.max(np.abs(G[:, 2, 2] - Gx[:, 2, 2]) / np.abs(Gx[:, 2, 2]))
            conv += [(key, "in-plane", m, err_in), (key, "out-of-plane", m, err_out)]
        print(f"{key}: max rel error of the defect-site G vs m: " +
              ", ".join(f"{m}: {e_i:.1e}/{e_o:.1e}" for (_, _, m, e_i), (_, _, _, e_o) in zip(conv[-2 * len(M_SCAN)::2], conv[-2 * len(M_SCAN) + 1::2])))
        G = cluster_green_cf(Dcc, B0, A, Bc, (w_up + 1j * ETA_D) ** 2)[:, own][:, :, own]
        ldos_up += [(key, w, a, b) for w, a, b in zip(w_up, *rho(w_up, G))]
        w_ex = (np.r_[np.arange(0.05, 3.3, 0.05), [1.44, 1.45, 1.455, 1.46, 3.14, 3.15, 3.155, 3.16, 3.17]] if key == "light"
                else np.r_[np.arange(0.01, 0.3, 0.01), np.arange(0.3, 2.0, 0.05)])
        w_ex = np.unique(np.round(w_ex, 4))
        ldos_ex += [(key, w, a, b) for w, a, b in zip(w_ex, *rho(w_ex, site_green(D, 3 * site + np.arange(3), w_ex, ETA_D)))]
        if key == "light":
            for lab, target, k in LOCAL:
                exact = np.sqrt(np.sort(sla.eigsh(D, k=k, sigma=target ** 2, which="LM", return_eigenvectors=False)))
                for m in M_LOCAL:
                    w2 = np.linalg.eigvalsh(assemble(Dcc, B0, A[:m], Bc[:m - 1]))
                    up = np.sqrt(np.clip(w2, 0, None))
                    near = up[np.argmin(np.abs(up - exact[0]))]
                    local.append((lab, float(exact[0]), m, float(near), abs(float(near - exact[0]))))
                print(f"light {lab}: exact {exact}; upfolded |error| by m: " +
                      ", ".join(f"{m}: {e:.1e}" for (_, _, m, _, e) in local[-len(M_LOCAL):]))
        print(f"{key} done ({time.time() - t0:.0f} s)")
        write_defect_csv(conv, local, ldos_up, ldos_ex)


def write_defect_csv(conv, local, ldos_up, ldos_ex):
    write_csv("honeycomb-convergence.csv", ["defect", "polarization", "chain_blocks_m", "max_rel_err_G00"],
              [(k, p, m, f"{e:.4e}") for k, p, m, e in conv])
    write_csv("honeycomb-local-modes.csv", ["mode", "exact_omega", "chain_blocks_m", "upfolded_omega", "abs_err"],
              [(l, f"{x:.10f}", m, f"{u:.10f}", f"{e:.3e}") for l, x, m, u, e in local])
    write_csv("honeycomb-ldos-upfolded.csv", ["defect", "omega", "rho_x", "rho_z"],
              [(k, f"{w:.5f}", f"{a:.6e}", f"{b:.6e}") for k, w, a, b in ldos_up])
    write_csv("honeycomb-ldos-exact.csv", ["defect", "omega", "rho_x", "rho_z"],
              [(k, f"{w:.5f}", f"{a:.6e}", f"{b:.6e}") for k, w, a, b in ldos_ex])


def read(name):
    with open(os.path.join(DATA, name)) as f:
        return list(csv.DictReader(f))


def defect_figures():
    up, ex = read("honeycomb-ldos-upfolded.csv"), read("honeycomb-ldos-exact.csv")
    gap = (1.3603, 1.5492)
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    for ax, (key, title, Md, fac) in zip(axes, DEFECTS):
        for pol, col, lab in (("x", COLOR["in"], "in-plane (x)"), ("z", COLOR["out"], "out-of-plane (z)")):
            w = np.array([float(r["omega"]) for r in up if r["defect"] == key])
            r_ = np.array([float(r["rho_" + pol]) for r in up if r["defect"] == key])
            ax.plot(w, np.clip(r_, 1e-4, None), color=col, lw=1.5, label=f"{lab}, upfolded")
            we = np.array([float(r["omega"]) for r in ex if r["defect"] == key])
            re_ = np.array([float(r["rho_" + pol]) for r in ex if r["defect"] == key])
            ax.plot(we, np.clip(re_, 1e-4, None), ls="none", marker="o", ms=5, mfc=col, mec="white", mew=1.2,
                    label=f"{lab}, exact supercell", zorder=5)
        ax.axvspan(*gap, color=AXIS, alpha=0.25, lw=0)
        ax.axvline(1.8621, color=AXIS, lw=0.75)
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 100)
        ax.set_xlim(0, 3.3 if key == "light" else 2.0)
        ax.set_xlabel("Frequency ω")
        ax.set_title(f"{title} (M′ = {Md:g}, springs × {fac:g})", loc="left", color=INK, fontsize=11)
    axes[0].set_ylabel("Local density of states at the defect")
    axes[0].text(0.5 * sum(gap), 40, "gap", ha="center", va="center", color=INK2, fontsize=9)
    axes[0].text(1.8621 + 0.04, 40, "host top", ha="left", va="center", color=INK2, fontsize=9)
    arrow = dict(arrowstyle="-", color=MUTED, lw=0.75, shrinkB=3)
    axes[0].annotate("out-of-plane mode\nin the gap, 1.454", xy=(1.42, 8.0), xytext=(1.25, 13), ha="right", va="center",
                     color=INK2, fontsize=9, arrowprops=arrow)
    axes[0].annotate("in-plane local\nmode, 3.155", xy=(3.12, 9.0), xytext=(2.9, 18), ha="right", va="center",
                     color=INK2, fontsize=9, arrowprops=arrow)
    axes[1].annotate("flexural resonance", xy=(0.075, 5.5), xytext=(0.32, 22), ha="left", va="center",
                     color=INK2, fontsize=9, arrowprops=arrow)
    axes[1].annotate("in-plane resonance", xy=(0.29, 3.6), xytext=(0.55, 7), ha="left", va="center",
                     color=INK2, fontsize=9, arrowprops=arrow)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=4, fontsize=9, labelcolor=INK2)
    save(fig, "honeycomb-ldos")

    conv, loc = read("honeycomb-convergence.csv"), read("honeycomb-local-modes.csv")
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.0, 3.8), constrained_layout=True)
    for key, title, _, _ in DEFECTS:
        for pol, col in (("in-plane", COLOR["in"]), ("out-of-plane", COLOR["out"])):
            rows = [r for r in conv if r["defect"] == key and r["polarization"] == pol]
            ax.plot([int(r["chain_blocks_m"]) for r in rows], [float(r["max_rel_err_G00"]) for r in rows], color=col, lw=1.5,
                    marker="o" if key == "light" else "s", ms=6, mec="white", mew=1.2, label=f"{title}, {pol}")
    ax.set_yscale("log")
    ax.set_xlabel("Chain blocks m")
    ax.set_ylabel("Max relative error of G at the defect")
    ax.set_title("Defect-site Green function against the supercell", loc="left", color=INK, fontsize=11)
    ax.legend(fontsize=9, labelcolor=INK2)
    for (lab, _, _), col in zip(LOCAL[:2], (COLOR["in"], COLOR["out"])):      # the deep bound states
        rows = [r for r in loc if r["mode"] == lab and int(r["chain_blocks_m"]) <= 12]
        ax2.plot([int(r["chain_blocks_m"]) for r in rows], [max(float(r["abs_err"]), 1e-16) for r in rows], color=col, lw=1.5,
                 marker="o", ms=6, mec="white", mew=1.2, label=f"{lab} ({float(rows[0]['exact_omega']):.4f})")
    ax2.set_yscale("log")
    ax2.set_xlabel("Chain blocks m")
    ax2.set_ylabel("|ω upfolded − ω exact|")
    ax2.set_title("Light defect: bound states", loc="left", color=INK, fontsize=11)
    ax2.legend(fontsize=9, labelcolor=INK2)
    save(fig, "honeycomb-convergence")


# ---------------------------------------------------------------- polar switch
L_POLAR, MMAX_P = 32, 200
M_SCAN_P = [5, 10, 20, 40, 60, 80, 120, 160, 200]


def polar_data():
    """Light and heavy defects on the non-polar and the polar host, 32 x 32 mesh: chain from the Bloch
    modes, exact Green functions from one dense diagonalization of each defect supercell; for the polar
    host also the upfolded result without the column scaling of B0."""
    t0 = time.time()
    L = L_POLAR
    site = hc.site_index(0, 0, 0, L)
    hop = hc.hop_distance(site, L)
    C = hc.coordinates(np.where(hop <= RADIUS)[0])
    own = [int(np.where(C == 3 * site + a)[0][0]) for a in range(3)]
    z_test = (TEST + 1j * ETA_D) ** 2
    w_up = np.linspace(0.005, 3.3, 1320)
    conv, ldos = [], []
    for host, polar in (("non-polar", None), ("polar", hc.POLAR)):
        B0, A, Bc, (A1, W) = hc.mode_chain(C, L, MMAX_P, polar=polar)
        D0, m0 = hc.supercell(L, polar=polar)
        D0 = D0.toarray() if sp.issparse(D0) else D0
        print(f"{host}: chain {len(A)} blocks, widths {A[0].shape[0]}..{A[-1].shape[0]}, first block "
              f"{np.abs(A1 - W.conj().T @ D0[np.ix_(C, C)] @ W).max():.1e} ({time.time() - t0:.0f} s)", flush=True)
        for key, _, Md, fac in DEFECTS:
            D, m = hc.supercell(L, defect=(0, (0, 0), Md, fac), polar=polar)
            D = D.toarray() if sp.issparse(D) else D
            lam, V = np.linalg.eigh(D)
            Vd = V[3 * site:3 * site + 3]

            def exact(z):
                return np.einsum("an,zn,bn->zab", Vd, 1.0 / (z[:, None] - lam[None, :]), Vd)
            Gx = exact(z_test)
            Dcc = D[np.ix_(C, C)]
            variants = [("scaled", hc.defect_coupling(B0, C, m0, m))]
            if polar is not None:
                variants.append(("naive", B0))
            for name, b0 in variants:
                for mm in M_SCAN_P:
                    G = cluster_green_cf(Dcc, b0, A[:mm], Bc[:mm - 1], z_test)[:, own][:, :, own]
                    e_in = max(np.max(np.abs(G[:, a, a] - Gx[:, a, a]) / np.abs(Gx[:, a, a])) for a in (0, 1))
                    e_out = np.max(np.abs(G[:, 2, 2] - Gx[:, 2, 2]) / np.abs(Gx[:, 2, 2]))
                    conv += [(host, key, name, "in-plane", mm, e_in), (host, key, name, "out-of-plane", mm, e_out)]
                print(f"  {host} {key} {name}: " + ", ".join(f"{c[4]}: {c[5]:.1e}" for c in conv[-2 * len(M_SCAN_P)::2]), flush=True)
            zz = (w_up + 1j * ETA_D) ** 2
            G = cluster_green_cf(Dcc, variants[0][1], A, Bc, zz)[:, own][:, :, own]
            ldos += [(host, key, "upfolded", w, a, b) for w, a, b in zip(w_up, *rho(w_up, G))]
            w_ex = w_up[::22]
            ldos += [(host, key, "exact", w, a, b) for w, a, b in zip(w_ex, *rho(w_ex, exact((w_ex + 1j * ETA_D) ** 2)))]
            write_csv("honeycomb-polar-convergence.csv", ["host", "defect", "coupling", "polarization", "chain_blocks_m", "max_rel_err_G00"],
                      [(h, k, n, p_, mm, f"{e:.4e}") for h, k, n, p_, mm, e in conv])
            write_csv("honeycomb-polar-ldos.csv", ["host", "defect", "source", "omega", "rho_x", "rho_z"],
                      [(h, k, src, f"{w:.5f}", f"{a:.6e}", f"{b:.6e}") for h, k, src, w, a, b in ldos])
            print(f"  {host} {key} done ({time.time() - t0:.0f} s)", flush=True)


def polar_figures():
    # dispersion with and without the long-range term
    s, kpath, ticks = hc.path(per_unit=80)
    fig, (ax, axz) = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.6, 1.0]))
    for polar, ls, lab in ((None, (0, (4, 3)), "non-polar"), (hc.POLAR, "-", "polar")):
        w = np.array([np.sqrt(np.clip(np.linalg.eigvalsh(hc.bloch(k, polar=polar)), 0, None)) for k in kpath])
        for n in range(6):
            ax.plot(s, w[:, n], color=INK2 if polar is None else COLOR["in"], lw=1.0 if polar is None else 1.5, ls=ls,
                    label=lab if n == 0 else None)
    ax.set_xlim(0, ticks[-1])
    ax.set_xticks(ticks, ["Γ", "M", "K", "Γ"])
    ax.grid(axis="x", visible=False)
    for t in ticks[1:-1]:
        ax.axvline(t, color=AXIS, lw=0.75)
    ax.set_ylim(0, 2.0)
    ax.set_ylabel("Frequency ω")
    ax.set_title("Dispersion with and without the polar term", loc="left", color=INK, fontsize=11)
    ax.legend(loc="lower center", ncol=2, fontsize=9.5, labelcolor=INK2)
    q = np.linspace(0, 0.6, 241)
    for direction, col, lab in ((np.array([1.0, 0.0]), COLOR["in"], "along Γ–K"), (hc.POINTS["M"] / np.linalg.norm(hc.POINTS["M"]), COLOR["out"], "along Γ–M")):
        w = np.array([np.sqrt(np.clip(np.linalg.eigvalsh(hc.bloch(qq * direction, polar=hc.POLAR)), 0, None)) for qq in q])
        axz.plot(q, w[:, 5] - w[:, 4], color=col, lw=1.5, label=lab)
    axz.plot(q, 0.302 * q, color=INK2, lw=0.75, ls=(0, (4, 3)), label="0.302 |q|")
    axz.set_xlim(0, 0.6)
    axz.set_ylim(0, None)
    axz.set_xlabel("|q|")
    axz.set_ylabel("ω_LO − ω_TO")
    axz.set_title("LO–TO splitting near Γ", loc="left", color=INK, fontsize=11)
    axz.legend(fontsize=9.5, labelcolor=INK2)
    save(fig, "honeycomb-polar-dispersion")

    rows = read("honeycomb-polar-convergence.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8), constrained_layout=True, sharey=True)
    for ax, (key, title, _, _) in zip(axes, DEFECTS):
        for host, coupling, col, mk, lab in (("non-polar", "scaled", INK2, "o", "non-polar host"),
                                             ("polar", "scaled", COLOR["in"], "o", "polar host"),
                                             ("polar", "naive", COLOR["out"], "s", "polar host, B₀ not scaled")):
            for pol, ls in (("in-plane", "-"),):
                rr = [r for r in rows if r["host"] == host and r["defect"] == key and r["coupling"] == coupling and r["polarization"] == pol]
                ax.plot([int(r["chain_blocks_m"]) for r in rr], [float(r["max_rel_err_G00"]) for r in rr], color=col, lw=1.5, ls=ls,
                        marker=mk, ms=6, mec="white", mew=1.2, label=lab)
        ax.set_yscale("log")
        ax.set_xlabel("Chain blocks m")
        ax.set_title(f"{title}: in-plane G at the defect", loc="left", color=INK, fontsize=11)
    axes[0].set_ylabel("Max relative error against the supercell")
    axes[0].legend(fontsize=9, labelcolor=INK2)
    save(fig, "honeycomb-polar-convergence")


# ---------------------------------------------------------------- spectral function at a finite concentration
SPEC_LQ, SPEC_ETA, SPEC_ND = 384, 0.01, 0.02       # host mesh, broadening, defects per cell (on A sites)
MVEC = hc.masses()
FREE = [(1.3603, 1.5492), (1.8621, np.inf)]          # frequencies without host modes: gap and above the top


def averaged_spectral(qs, omegas, eta, T, offs, subl, nd, polar=None):
    """Configuration-averaged A(q, w) = -(2w/pi) Im Tr[M^1/2 <g> M^1/2] (host masses), split into in-plane and
    out-of-plane, with Sigma(q) from `hc.self_energy`. Dyson equation where the host has modes; first order,
    g0 + g0 Sigma g0, where it has none (FREE), which keeps bound states of a single defect at their
    frequencies (checked against random configurations, examples/honeycomb_spectral_check.py)."""
    om = np.asarray(omegas)
    z = (om + 1j * eta) ** 2
    free = np.zeros(len(om), bool)
    for lo, hi in FREE:
        free |= (om > lo) & (om < hi)
    Wm = np.sqrt(np.outer(MVEC, MVEC))
    A_in, A_out = np.empty((len(qs), len(om))), np.empty((len(qs), len(om)))
    for k, q in enumerate(qs):
        Sig = hc.self_energy(q, T, offs, subl, nd)
        Phiq = hc.bloch(q, gauge="cell", polar=polar) * Wm
        g0 = np.linalg.inv(z[:, None, None] * np.diag(MVEC)[None] - Phiq[None])
        g = np.where(free[:, None, None], g0 + g0 @ Sig @ g0, np.linalg.inv(np.linalg.inv(g0) - Sig))
        tr = -(2 * om[:, None] / np.pi) * np.einsum("wii->wi", g * Wm).imag
        A_in[k], A_out[k] = tr[:, hc.IN_PLANE].sum(1), tr[:, hc.OUT_OF_PLANE].sum(1)
    return A_in, A_out


def spectral_data():
    t0 = time.time()
    s, kpath, ticks = hc.path(per_unit=60)
    om = np.arange(0.002, 3.4, 0.004)
    omz, kz = np.arange(0.001, 0.45, 0.001), np.linspace(0.0, 0.8, 161)
    qz = np.outer(kz, hc.POINTS["M"] / np.linalg.norm(hc.POINTS["M"]))      # along Γ-M
    out, cuts = {}, []
    offs, subl, _, _ = hc.perturbation((0, (0, 0), 0.3, 1.3))
    g_main = hc.host_green_support(om, SPEC_ETA, SPEC_LQ, offs, subl)
    g_zoom = hc.host_green_support(omz, SPEC_ETA, SPEC_LQ, offs, subl)
    print(f"host Green functions on the support ({time.time() - t0:.0f} s)", flush=True)
    for key, _, Md, fac in DEFECTS:
        offs, subl, dPhi, dM = hc.perturbation((0, (0, 0), Md, fac))
        T = hc.tmatrix_u(om, SPEC_ETA, g_main, dPhi, dM)
        out[key] = sum(averaged_spectral(kpath, om, SPEC_ETA, T, offs, subl, SPEC_ND))
        for lab, q in (("Γ+0.05M", 0.05 * hc.POINTS["M"]), ("M", hc.POINTS["M"]), ("K", hc.POINTS["K"])):
            a_in, a_out = averaged_spectral(q[None], om, SPEC_ETA, T, offs, subl, SPEC_ND)
            cuts += [(key, lab, w, x, y) for w, x, y in zip(om, a_in[0], a_out[0])]
        if key == "heavy":
            Tz = hc.tmatrix_u(omz, SPEC_ETA, g_zoom, dPhi, dM)
            out["zoom"] = sum(averaged_spectral(qz, omz, SPEC_ETA, Tz, offs, subl, SPEC_ND))
        print(f"{key} done ({time.time() - t0:.0f} s)", flush=True)
    cache = os.environ.get("HONEYCOMB_SPECTRAL_CACHE")       # optional .npz outside the repository, for --spectral-plot
    if cache:
        np.savez(cache, s=s, ticks=ticks, om=om, omz=omz, kz=kz, **out)
    write_csv("honeycomb-spectral-cuts.csv", ["defect", "q", "omega", "A_in_plane", "A_out_of_plane"],
              [(k, l, f"{w:.4f}", f"{a:.6e}", f"{b:.6e}") for k, l, w, a, b in cuts])
    return s, ticks, om, omz, kz, out


def spectral_figure(s, ticks, om, omz, kz, out):
    fig = plt.figure(figsize=(10.0, 4.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 4, width_ratios=[1.0, 1.0, 0.68, 0.045])
    axa = fig.add_subplot(gs[0])
    axb = fig.add_subplot(gs[1], sharey=axa)
    axc, cax = fig.add_subplot(gs[2]), fig.add_subplot(gs[3])
    levels, colors = [0.1, 0.3, 1.0, 3.0, 10.0, 1e6], ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
    for ax, key, title in ((axa, "light", "Light defect"), (axb, "heavy", "Heavy defect")):
        cs = ax.contourf(s[::2], om, np.clip(out[key][::2].T, 1e-12, None), levels=levels, colors=colors)
        ax.grid(False)
        ax.set_xticks(ticks, ["Γ", "M", "K", "Γ"])
        for t in ticks[1:-1]:
            ax.axvline(t, color=AXIS, lw=0.75)
        ax.axhline(1.8621, color=AXIS, lw=0.75)
        ax.set_xlim(0, ticks[-1])
        ax.set_ylim(0, 3.4)
        ax.set_title(f"{title}, $n_d$ = {SPEC_ND} per cell", loc="left", color=INK, fontsize=11)
    axb.tick_params(labelleft=False)
    axa.set_ylabel("Frequency ω")
    arrow = dict(arrowstyle="-", color=MUTED, lw=0.75, shrinkB=2)
    axa.annotate("in-plane local modes, 3.155", xy=(ticks[2] - 0.5, 3.155 - 0.03), xytext=(ticks[2] - 1.3, 2.55),
                 ha="center", va="center", color=INK2, fontsize=9, arrowprops=arrow)
    axa.annotate("out-of-plane modes\nin the gap, 1.454", xy=(ticks[1] * 0.45, 1.454), xytext=(ticks[1] * 0.45, 2.25),
                 ha="center", va="center", color=INK2, fontsize=9, arrowprops=arrow)
    axb.add_patch(matplotlib.patches.Rectangle((0.0, 0.0), kz[-1], omz[-1], fill=False, ec=INK2, lw=0.8))
    axb.text(kz[-1] + 0.08, 0.04, "zoom", ha="left", va="bottom", color=INK2, fontsize=9)
    axc.contourf(kz, omz, np.clip(out["zoom"].T, 1e-12, None), levels=levels, colors=colors)
    axc.grid(False)
    qd = np.outer(kz, hc.POINTS["M"] / np.linalg.norm(hc.POINTS["M"]))
    wb = np.array([np.sqrt(np.clip(np.linalg.eigvalsh(hc.bloch(q)), 0, None)) for q in qd])
    for n in range(3):
        axc.plot(kz, wb[:, n], color=INK2, lw=0.8, ls=(0, (4, 3)))
    axc.axhline(0.053, color=INK2, lw=0.75)
    axc.axhline(0.265, color=INK2, lw=0.75)
    box = dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.85)
    axc.text(kz[-1] - 0.02, 0.053 + 0.008, "flexural resonance, 0.053", ha="right", va="bottom", color=INK2, fontsize=8.5, bbox=box)
    axc.text(kz[-1] - 0.02, 0.265 + 0.008, "in-plane resonance, 0.265", ha="right", va="bottom", color=INK2, fontsize=8.5, bbox=box)
    axc.text(0.02, omz[-1] - 0.01, "dashed: host ZA, TA, LA", ha="left", va="top", color=INK2, fontsize=8.5, bbox=box)
    axc.set_xlim(0, kz[-1])
    axc.set_ylim(0, omz[-1])
    axc.set_xlabel("q along Γ–M")
    axc.set_title("Heavy defect near Γ", loc="left", color=INK, fontsize=11)
    cb = fig.colorbar(cs, cax=cax, ticks=levels[:-1])
    cb.ax.set_yticklabels(["0.1", "0.3", "1", "3", "10"])
    cb.set_label("A(q, ω)", color=INK2)
    cb.outline.set_visible(False)
    save(fig, "honeycomb-spectral")


def main():
    ap = argparse.ArgumentParser()
    for flag in ("host", "defects", "plot", "polar", "polar_plot", "spectral", "spectral_plot"):
        ap.add_argument("--" + flag.replace("_", "-"), action="store_true")
    args = ap.parse_args()
    if args.polar:
        polar_data()
    if args.polar or args.polar_plot:
        polar_figures()
    if args.spectral:
        spectral_figure(*spectral_data())
    if args.spectral_plot:
        z = np.load(os.environ["HONEYCOMB_SPECTRAL_CACHE"])
        spectral_figure(z["s"], list(z["ticks"]), z["om"], z["omz"], z["kz"], {k: z[k] for k in ("light", "heavy", "zoom")})
    run_all = not (args.host or args.defects or args.plot or args.polar or args.polar_plot or args.spectral or args.spectral_plot)
    if args.host or run_all:
        host_figure()
    if args.defects or run_all:
        defect_data()
    if args.plot or args.defects or run_all:
        defect_figures()


if __name__ == "__main__":
    main()
