#!/usr/bin/env python3
"""Figures and data for the toy-validation page (run: npm run figures, about one minute).

Writes site/results/toy-local-dos.{svg,pdf,png}, site/results/toy-convergence.{svg,pdf,png}
and the plotted numbers as CSV under site/data/. Same lattice, defects, cluster and
chain as examples/square_lattice_toy.py.

Colour: the chain length m is ordered, so the upfolded curves use one blue ramp
(light = short chain); the exact reference is neutral ink. The two defects in the
convergence figure are categories (palette slots 1 and 2). Both palettes were checked
with the data-viz validator against the white page surface.
"""
import csv
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from phdef import block_lanczos, assemble, modes, cluster_green_modes  # noqa: E402
from phdef.models import square_lattice_dynmat, manhattan_partition  # noqa: E402

L, RC, MMAX, ETA = 121, 2, 160, 0.05
DEFECTS = [("light", "Light defect (M′ = 0.3, K′ = 1.3)", 0.3, 1.3),
           ("heavy", "Heavy defect (M′ = 8, K′ = 0.6)", 8.0, 0.6)]
M_CURVES = [10, 40, 160]
RAMP = {10: "#86b6ef", 40: "#2a78d6", 160: "#104281"}          # ordinal blue, validated
CATEGORICAL = {"light": "#2a78d6", "heavy": "#eb6834"}          # slots 1-2, validated
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
TEST_OMEGAS = np.array([0.15, 0.25, 0.35, 0.5, 0.8, 1.2, 1.8, 2.5])
M_SCAN = [5, 10, 20, 40, 60, 80, 100, 120, 140, 160]
M_LOCAL = list(range(1, 11))
FLOOR = 1e-16                                                   # plotted floor for an exact zero

plt.rcParams.update({
    "svg.hashsalt": "phonon-defect", "svg.fonttype": "path", "pdf.fonttype": 42,
    "font.family": "DejaVu Sans", "font.size": 10.5,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.75, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.75, "grid.linestyle": "-",
    "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white",
})


def rho_from_green(omega, G):
    return -2.0 * omega / np.pi * G.imag


def exact_green_00(D, c0, omegas):
    e0 = np.zeros(D.shape[0], dtype=complex)
    e0[c0] = 1.0
    Dcsc = D.tocsc()
    I = sp.identity(D.shape[0], format="csc")
    return np.array([sla.spsolve(((w + 1j * ETA) ** 2) * I - Dcsc, e0)[c0] for w in omegas])


def save(fig, name):
    out = os.path.join(ROOT, "site", "results")
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, name + ".svg"), metadata={"Date": None})
    fig.savefig(os.path.join(out, name + ".pdf"), metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(os.path.join(out, name + ".png"), dpi=200)
    plt.close(fig)


def write_csv(name, header, rows):
    out = os.path.join(ROOT, "site", "data")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def main():
    C, Bk = manhattan_partition(L, RC)
    host, c0 = square_lattice_dynmat(L)
    i0 = int(np.where(C == c0)[0][0])
    B0, A_all, B_all = block_lanczos(host[Bk][:, Bk].tocsr(), host[Bk][:, C].toarray(), MMAX)

    curves, exact_pts, conv = {}, {}, {}
    for key, _, Md, Kd in DEFECTS:
        D, _ = square_lattice_dynmat(L, Md, Kd)
        D_CC = D[C][:, C].toarray()
        wmax = 4.6 if key == "light" else 3.0
        w_dense = np.linspace(0.02, wmax, 900)
        if key == "light":
            w_exact = np.concatenate([np.arange(0.1, 4.25, 0.1), np.arange(4.26, 4.451, 0.02)])
        else:
            w_exact = np.arange(0.1, 3.0, 0.1)
        exact_pts[key] = (w_exact, rho_from_green(w_exact, exact_green_00(D, c0, w_exact)))
        curves[key] = {}
        top_exact = np.sqrt(sla.eigsh(D, k=1, which="LA", return_eigenvectors=False)[0])
        G_test_exact = exact_green_00(D, c0, TEST_OMEGAS)
        conv[key] = []
        for m in M_SCAN:
            w2, V = modes(assemble(D_CC, B0, A_all[:m], B_all[:m - 1]))
            G_test = cluster_green_modes(w2, V, len(C), (TEST_OMEGAS + 1j * ETA) ** 2)[:, i0, i0]
            err = float(np.max(np.abs(G_test - G_test_exact) / np.abs(G_test_exact)))
            top_err = abs(float(np.sqrt(w2[-1])) - float(top_exact))
            conv[key].append((m, err, top_err))
            if m in M_CURVES:
                G = cluster_green_modes(w2, V, len(C), (w_dense + 1j * ETA) ** 2)[:, i0, i0]
                curves[key][m] = (w_dense, rho_from_green(w_dense, G))
        print(f"{key}: exact top mode {top_exact:.6f}; " +
              ", ".join(f"m={m}: err {e:.1e}" for m, e, _ in conv[key]))
        if key == "light":
            local_mode = []
            for m in M_LOCAL:
                w2, _ = modes(assemble(D_CC, B0, A_all[:m], B_all[:m - 1]))
                local_mode.append((m, abs(float(np.sqrt(w2[-1])) - float(top_exact))))
            print("light local mode |error| by m:", ", ".join(f"{m}: {e:.1e}" for m, e in local_mode))

    # ---- figure 1: local density of states at the defect ----
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    for ax, (key, title, _, _) in zip(axes, DEFECTS):
        for m in M_CURVES:
            w, r = curves[key][m]
            ax.plot(w, np.clip(r, 1e-4, None), color=RAMP[m], lw=1.5, solid_capstyle="round",
                    label=f"Upfolded, m = {m} chain blocks")
        we, re = exact_pts[key]
        ax.plot(we, re, ls="none", marker="o", ms=6, mfc=INK2, mec="white", mew=1.5,
                label="Exact (full lattice)", zorder=5)
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 30)
        ax.set_xlabel("Frequency ω")
        ax.set_title(title, loc="left", color=INK, fontsize=11)
        ax.axvline(2 * np.sqrt(2), color=AXIS, lw=0.75)
        if key == "light":     # the light panel has its local-mode label up top; the curve is high here
            ax.text(2 * np.sqrt(2) - 0.05, 1.3e-3, "host band top", rotation=90, ha="right", va="bottom",
                    color=MUTED, fontsize=9)
        else:                  # the heavy panel's curve runs low near the band top
            ax.text(2 * np.sqrt(2) - 0.05, 22, "host band top", rotation=90, ha="right", va="top",
                    color=MUTED, fontsize=9)
    axes[0].set_ylabel("Local density of states ρ₀₀(ω)")
    axes[0].annotate("local mode ω = 4.3488,\nexact for every m", xy=(4.29, 6.0), xytext=(4.05, 2.0),
                     color=INK2, fontsize=9, ha="right", va="center", zorder=6,
                     bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"),
                     arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.75, shrinkA=2, shrinkB=4))
    axes[1].annotate("in-band resonance", xy=(0.5, 2.3), xytext=(1.0, 6.0),
                     color=INK2, fontsize=9, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.75))
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=4, fontsize=9.5, labelcolor=INK2)
    save(fig, "toy-local-dos")

    # ---- figure 2: convergence with the chain length ----
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8), constrained_layout=True)
    ax = axes[0]
    for key, title, _, _ in DEFECTS:
        ms = [c[0] for c in conv[key]]
        ax.plot(ms, [c[1] for c in conv[key]], color=CATEGORICAL[key], lw=1.5, marker="o", ms=6,
                mec="white", mew=1.5, label=title)
    ax.set_yscale("log")
    ax.set_xlabel("Chain blocks m")
    ax.set_ylabel("Max relative error of G₀₀(ω + iη)")
    ax.set_title("In-band response: error falls as the chain grows", loc="left", color=INK, fontsize=11)
    ax.legend(fontsize=9.5, labelcolor=INK2)
    ax = axes[1]
    ax.plot([m for m, _ in local_mode], [max(e, FLOOR) for _, e in local_mode], color=CATEGORICAL["light"],
            lw=1.5, marker="o", ms=6, mec="white", mew=1.5)
    ax.set_yscale("log")
    ax.set_ylim(3e-17, 1)
    ax.set_xticks(M_LOCAL)
    ax.axhline(4.35 * np.finfo(float).eps, color=AXIS, lw=0.75)
    ax.text(10, 4.35 * np.finfo(float).eps * 2.5, "double-precision floor", ha="right", va="bottom",
            color=MUTED, fontsize=9)
    ax.set_xlabel("Chain blocks m")
    ax.set_ylabel("Local-mode error |ω upfolded − ω exact|")
    ax.set_title("Light defect: the local mode converges in a few blocks", loc="left", color=INK, fontsize=11)
    save(fig, "toy-convergence")
    write_csv("toy-local-mode.csv", ["chain_blocks_m", "local_mode_abs_err"],
              [(m, f"{e:.4e}") for m, e in local_mode])

    write_csv("toy-local-dos-upfolded.csv", ["defect", "chain_blocks_m", "omega", "rho_00"],
              [(k, m, f"{w:.6f}", f"{r:.8e}") for k in curves for m in curves[k]
               for w, r in zip(*curves[k][m])])
    write_csv("toy-local-dos-exact.csv", ["defect", "omega", "rho_00"],
              [(k, f"{w:.6f}", f"{r:.8e}") for k in exact_pts for w, r in zip(*exact_pts[k])])
    write_csv("toy-convergence.csv", ["defect", "chain_blocks_m", "max_rel_err_G00", "top_mode_abs_err"],
              [(k, m, f"{e:.4e}", f"{t:.4e}") for k in conv for m, e, t in conv[k]])
    print("wrote site/results/toy-local-dos.*, site/results/toy-convergence.*, site/data/toy-*.csv")


if __name__ == "__main__":
    main()
