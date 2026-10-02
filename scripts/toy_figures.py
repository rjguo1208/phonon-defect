#!/usr/bin/env python3
"""Figures and data for the toy-validation page (run: npm run figures, about one minute).

Writes site/results/toy-local-dos.{svg,pdf,png}, site/results/toy-convergence.{svg,pdf,png},
site/results/toy-dispersion.{svg,pdf,png} (host dispersion beside the local densities of
states), site/results/toy-spectral.{svg,pdf,png} (A(k, w) at a defect concentration n_d in the
dilute limit, Sigma = n_d T, with the T-matrix in site/data/toy-tmatrix.csv) and the plotted
numbers as CSV under site/data/. Same lattice, defects, cluster
and chain as examples/square_lattice_toy.py.

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


PATH = [("Γ", 0.0, 0.0), ("X", np.pi, 0.0), ("M", np.pi, np.pi), ("Γ", 0.0, 0.0)]


def band_path(points_per_pi=160):
    """Host dispersion omega(k) = [2(2 - cos kx - cos ky)]^(1/2) along Γ-X-M-Γ (K = M = 1)."""
    ks, ss, ticks = [], [], [0.0]
    for n, ((_, x0, y0), (_, x1, y1)) in enumerate(zip(PATH[:-1], PATH[1:])):
        length = np.hypot(x1 - x0, y1 - y0)
        t = np.linspace(0.0, 1.0, int(round(points_per_pi * length / np.pi)) + 1)
        if n < len(PATH) - 2:
            t = t[:-1]                                      # the next segment starts at this corner
        ks.append(np.c_[x0 + (x1 - x0) * t, y0 + (y1 - y0) * t])
        ss.append(ticks[-1] + length * t)
        ticks.append(ticks[-1] + length)
    k = np.vstack(ks)
    return np.concatenate(ss), k, np.sqrt(2.0 * (2.0 - np.cos(k[:, 0]) - np.cos(k[:, 1]))), ticks


def dispersion_figure(s, w_band, ticks, w_dos, rho, w_local):
    """Host dispersion (left) beside the local densities of states at the centre site (right)."""
    fig, (ax, axd) = plt.subplots(1, 2, figsize=(10.0, 4.4), sharey=True, constrained_layout=True,
                                  gridspec_kw=dict(width_ratios=[1.85, 1.0]))
    top = 2 * np.sqrt(2)
    ax.plot(s, w_band, color=INK, lw=1.5, solid_capstyle="round", solid_joinstyle="round")
    ax.set_xlim(0, ticks[-1])
    ax.set_xticks(ticks, [p[0] for p in PATH])
    ax.grid(axis="x", visible=False)
    for t in ticks[1:-1]:
        ax.axvline(t, color=AXIS, lw=0.75)
    ax.set_ylim(0, 4.95)
    ax.set_ylabel("Frequency ω")
    ax.set_title("Host phonon dispersion", loc="left", color=INK, fontsize=11)
    for a in (ax, axd):
        a.axhline(top, color=AXIS, lw=0.75)
    ax.text(ticks[-1] - 0.12, top + 0.07, "host band top 2√2 ≈ 2.83 (at M)", ha="right", va="bottom",
            color=MUTED, fontsize=9)
    ax.axhline(w_local, color=CATEGORICAL["light"], lw=1.5)
    ax.text(ticks[-1] - 0.12, w_local + 0.08, f"light defect: local mode ω = {w_local:.3f}, no dispersion",
            ha="right", va="bottom", color=INK2, fontsize=9)
    ax.annotate("van Hove saddle at X", xy=(ticks[1], 2.0), xytext=(ticks[1] + 0.35, 1.3), color=INK2,
                fontsize=9, ha="left", va="center", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.75, shrinkB=2))
    ins = ax.inset_axes([0.03, 0.585, 0.17, 0.285])      # Brillouin zone and the path
    ins.set_facecolor("none")
    sq = np.pi * np.array([[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, -1]])
    ins.plot(sq[:, 0], sq[:, 1], color=AXIS, lw=0.75)
    ins.plot([0, np.pi, np.pi, 0], [0, 0, np.pi, 0], color=INK, lw=1.2, solid_joinstyle="round")
    for lab, x, y, ha, va in (("Γ", -0.3, -0.3, "right", "top"), ("X", np.pi + 0.35, 0.0, "left", "center"),
                              ("M", np.pi + 0.35, np.pi, "left", "center")):
        ins.text(x, y, lab, ha=ha, va=va, color=INK2, fontsize=9)
    ins.set_xlim(-4.0, 5.2)
    ins.set_ylim(-4.0, 4.0)
    ins.set_aspect("equal")
    ins.axis("off")

    for key, label, color in (("host", "no defect (host)", INK2), ("light", "light defect", CATEGORICAL["light"]),
                              ("heavy", "heavy defect", CATEGORICAL["heavy"])):
        axd.plot(np.clip(rho[key], 1e-4, None), w_dos, color=color, lw=1.5, label=label, solid_capstyle="round")
    axd.set_xscale("log")
    axd.set_xlim(1e-3, 30)
    axd.set_xlabel("ρ₀₀(ω) at the centre site")
    axd.set_title("Local density of states", loc="left", color=INK, fontsize=11)
    axd.tick_params(labelleft=False)
    j = int(np.argmax(rho["light"]))
    axd.annotate("local mode", xy=(rho["light"][j], w_dos[j]), xytext=(1.6, w_dos[j] + 0.33), color=INK2,
                 fontsize=9, ha="center", va="bottom", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.75, shrinkB=2))
    jh = int(np.argmax(np.where(w_dos < 1.5, rho["heavy"], 0.0)))
    axd.annotate("resonance", xy=(rho["heavy"][jh], w_dos[jh]), xytext=(2.6, w_dos[jh] + 0.55), color=INK2,
                 fontsize=9, ha="center", va="bottom", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.75, shrinkB=2))
    axd.legend(loc="upper right", bbox_to_anchor=(1.0, 0.79), fontsize=9.5, labelcolor=INK2)
    save(fig, "toy-dispersion")


# ---- spectral function at a finite defect concentration (dilute limit, Sigma = n_d T) ----
SUPPORT = np.array([(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)])   # sites where V(z) acts
SHELLS = [(0, 0), (1, 0), (1, 1), (2, 0)]                          # distinct |r_a - r_b| on SUPPORT
ND, ETA_A = 0.02, 0.01                                             # concentration and broadening of A(k, w)
SEQ = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # sequential blue, steps 100/250/400/550/700
A_LEVELS = [0.1, 0.3, 1.0, 3.0, 10.0, 1e6]


def host_green(omegas, eta, nky=32768, chunk=64):
    """Infinite-lattice host Green function g0(r; z) = (1/N) sum_k e^{ik.r} / (z - w_k^2), z = (w + i eta)^2,
    for r in SHELLS: the k_x integral is done exactly, the k_y integral by the trapezoid rule on [0, pi]."""
    ky = np.linspace(0.0, np.pi, nky + 1)
    wts = np.full(nky + 1, 1.0 / nky)
    wts[[0, -1]] *= 0.5
    om = np.atleast_1d(np.asarray(omegas, float))
    out = np.empty((len(om), len(SHELLS)), complex)
    for i in range(0, len(om), chunk):
        w = (om[i:i + chunk, None] + 1j * eta) ** 2 - 4.0 + 2.0 * np.cos(ky)[None, :]
        r = (-w + np.sqrt(w * w - 4.0 + 0j)) / 2.0
        r = np.where(np.abs(r) < 1.0, r, 1.0 / r)              # root of x^2 + w x + 1 inside the unit circle
        base = 1.0 / (r - 1.0 / r)
        for j, (m, n) in enumerate(SHELLS):
            out[i:i + chunk, j] = (r ** m * base * np.cos(n * ky)[None, :]) @ wts
    return out


def defect_tmatrix(omegas, eta, Md, Kd, g):
    """Single-defect T-matrix on SUPPORT in displacement coordinates (host masses 1).

    V(z) = dPhi - z dM: the mass change at the defect and the spring change on its four bonds;
    T(z) = [1 - V(z) g0(z)]^-1 V(z). Returns T and the 5 x 5 host block g0 on SUPPORT."""
    z = (np.atleast_1d(np.asarray(omegas, float)) + 1j * eta) ** 2
    dK, dM = Kd - 1.0, Md - 1.0
    V = np.zeros((len(z), 5, 5), complex)
    V[:, 0, 0] = 4 * dK - z * dM
    V[:, 0, 1:] = -dK
    V[:, 1:, 0] = -dK
    V[:, np.arange(1, 5), np.arange(1, 5)] = dK
    idx = np.array([[SHELLS.index(tuple(sorted(np.abs(a - b), reverse=True))) for b in SUPPORT] for a in SUPPORT])
    G0 = g[:, idx]
    return np.linalg.solve(np.eye(5) - V @ G0, V), G0


def spectral_function(kpts, omegas, eta, T, nd, w_cut=np.inf):
    """Configuration-averaged A(k, w) = -(2w/pi) Im <g(k, z)> to first order in the concentration n_d.

    t_k(z) = sum_ab e^{-ik.r_a} T_ab(z) e^{ik.r_b}. Below w_cut: Dyson equation with Sigma = n_d t_k
    (the host branch needs the resummation). Above w_cut (the gap below a local mode): the first-order
    average g0 + n_d g0 t_k g0, which keeps the local mode at its isolated frequency; resumming n_d t_k
    near the pole of T shifts it by an artificial level repulsion that random positions do not produce."""
    om = np.asarray(omegas, float)
    phase = np.exp(1j * kpts @ SUPPORT.T)
    t = np.einsum("ka,wab,kb->kw", phase.conj(), T, phase)
    wk2 = 2 * (2 - np.cos(kpts[:, 0]) - np.cos(kpts[:, 1]))
    g0 = 1.0 / ((om[None, :] + 1j * eta) ** 2 - wk2[:, None])
    g = np.where(om[None, :] > w_cut, g0 + nd * g0 ** 2 * t, 1.0 / (1.0 / g0 - nd * t))
    return -(2 * om[None, :] / np.pi) * g.imag


def spectral_figure(s, ticks, om, A, kz, omz, Az, w_isolated, w_res):
    """A(k, w) along the path for both defects, and the heavy defect near Gamma along Gamma-X."""
    fig = plt.figure(figsize=(10.0, 4.5), constrained_layout=True)
    gs = fig.add_gridspec(1, 4, width_ratios=[1.0, 1.0, 0.68, 0.045])
    axa = fig.add_subplot(gs[0])
    axb = fig.add_subplot(gs[1], sharey=axa)
    axc, cax = fig.add_subplot(gs[2]), fig.add_subplot(gs[3])
    for ax, key, title in ((axa, "light", "Light defect"), (axb, "heavy", "Heavy defect")):
        cs = ax.contourf(s, om, np.clip(A[key].T, 1e-12, None), levels=A_LEVELS, colors=SEQ)
        ax.grid(False)
        ax.set_xticks(ticks, [p[0] for p in PATH])
        for t in ticks[1:-1]:
            ax.axvline(t, color=AXIS, lw=0.75)
        ax.axhline(2 * np.sqrt(2), color=AXIS, lw=0.75)
        ax.set_xlim(0, ticks[-1])
        ax.set_ylim(0, 4.9)
        ax.set_title(f"{title}, $n_d$ = {ND}", loc="left", color=INK, fontsize=11)
    axb.tick_params(labelleft=False)
    axa.set_ylabel("Frequency ω")
    axa.annotate(f"local modes: flat band at {w_isolated:.3f}", xy=(ticks[2] - 0.6, w_isolated - 0.03), xytext=(ticks[2] - 1.4, 3.6),
                 color=INK2, fontsize=9, ha="center", va="center",
                 arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.75, shrinkB=2))
    axb.add_patch(matplotlib.patches.Rectangle((0.0, 0.0), kz[-1], omz[-1], fill=False, ec=INK2, lw=0.8))
    axb.text(kz[-1] + 0.08, 0.04, "zoom", ha="left", va="bottom", color=INK2, fontsize=9)

    axc.contourf(kz, omz, np.clip(Az.T, 1e-12, None), levels=A_LEVELS, colors=SEQ)
    axc.grid(False)
    axc.plot(kz, 2 * np.sin(kz / 2), color=INK2, lw=0.8, ls=(0, (4, 3)))
    axc.text(kz[-1] - 0.03, 0.03, "dashed: host branch", ha="right", va="bottom", color=INK2, fontsize=9)
    axc.axhline(w_res, color=INK2, lw=0.75)
    axc.text(kz[-1] - 0.03, w_res + 0.02, f"one defect: {w_res:.2f}", ha="right", va="bottom", color=INK2, fontsize=9)
    axc.set_xlim(0, kz[-1])
    axc.set_ylim(0, omz[-1])
    axc.set_xlabel("k along Γ–X")
    axc.set_title("Heavy defect near Γ", loc="left", color=INK, fontsize=11)
    cb = fig.colorbar(cs, cax=cax, ticks=A_LEVELS[:-1])
    cb.ax.set_yticklabels(["0.1", "0.3", "1", "3", "10"])
    cb.set_label("A(k, ω)", color=INK2)
    cb.outline.set_visible(False)
    with plt.rc_context({"path.simplify_threshold": 0.5}):
        save(fig, "toy-spectral")


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

    curves, exact_pts, conv, modes_max = {}, {}, {}, {}
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
            if m == MMAX:
                modes_max[key] = (w2, V)
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
    modes_max["host"] = modes(assemble(host[C][:, C].toarray(), B0, A_all[:MMAX], B_all[:MMAX - 1]))

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

    # ---- figure 3: host dispersion beside the local densities of states (m = MMAX) ----
    s, kpath, w_band, ticks = band_path()
    w_dos = np.linspace(0.005, 4.95, 1980)
    z = (w_dos + 1j * ETA) ** 2
    rho = {}
    for key, (w2, V) in modes_max.items():
        rho[key] = rho_from_green(w_dos, (np.abs(V[i0]) ** 2 / (z[:, None] - w2[None, :])).sum(axis=1))
    w_local = float(np.sqrt(modes_max["light"][0][-1]))
    dispersion_figure(s, w_band, ticks, w_dos, rho, w_local)
    write_csv("toy-dispersion.csv", ["path_coordinate", "kx", "ky", "omega"],
              [(f"{a:.6f}", f"{kx:.6f}", f"{ky:.6f}", f"{w:.6f}") for a, (kx, ky), w in zip(s, kpath, w_band)])
    write_csv("toy-dos.csv", ["omega", "rho_00_no_defect", "rho_00_light", "rho_00_heavy"],
              [(f"{w:.6f}", f"{a:.8e}", f"{b:.8e}", f"{c:.8e}")
               for w, a, b, c in zip(w_dos, rho["host"], rho["light"], rho["heavy"])])
    jh = int(np.argmax(np.where(w_dos < 1.5, rho["heavy"], 0.0)))
    print(f"dispersion figure: local mode {w_local:.6f}; heavy resonance peak at {w_dos[jh]:.4f}; "
          f"host LDOS maximum at {w_dos[int(np.argmax(rho['host']))]:.4f}")

    # ---- figure 4: A(k, w) at a finite defect concentration, Sigma = n_d T ----
    om, kA = np.arange(0.002, 4.9, 0.005), slice(None, None, 2)       # every second k of the path
    omz, kz = np.arange(0.002, 1.2, 0.002), np.linspace(0.0, 1.3, 261)
    gA, gz = host_green(om, ETA_A), host_green(omz, ETA_A)
    A, T_out = {}, {}
    w_cut = 0.5 * (2 * np.sqrt(2) + w_local)                          # middle of the gap below the local mode
    for key, _, Md, Kd in DEFECTS:
        T_out[key], _ = defect_tmatrix(om, ETA_A, Md, Kd, gA)
        A[key] = spectral_function(kpath[kA], om, ETA_A, T_out[key], ND, w_cut)
        sums = [float(np.trapezoid(A[key][int(np.argmin(np.abs(s[kA] - x)))], om)) for x in (s[2], ticks[1], ticks[2])]
        print(f"A(k, w) {key}: weight near Γ / at X / at M {sums[0]:.4f} / {sums[1]:.4f} / {sums[2]:.4f}; "
              f"1 + n_d (1/M' - 1) = {1 + ND * (1 / Md - 1):.4f}")
    Md_h, Kd_h = [(d[2], d[3]) for d in DEFECTS if d[0] == "heavy"][0]
    Tz, G0z = defect_tmatrix(omz, ETA_A, Md_h, Kd_h, gz)
    Az = spectral_function(np.c_[kz, np.zeros_like(kz)], omz, ETA_A, Tz, ND)
    g00 = G0z[:, 0, 0] + np.einsum("wa,wab,wb->w", G0z[:, 0, :], Tz, G0z[:, :, 0])
    w_res = float(omz[np.argmax(-(2 * omz / np.pi) * (Md_h * g00).imag)])
    spectral_figure(s[kA], ticks, om, A, kz, omz, Az, w_local, w_res)
    comps = ((0, 0), (0, 1), (1, 1), (1, 2), (1, 3))
    write_csv("toy-tmatrix.csv", ["defect", "omega"] + [f"{p}_T{a}{b}" for a, b in comps for p in ("re", "im")],
              [(key, f"{w:.4f}", *[f"{v:.8e}" for a, b in comps for v in (T_out[key][i, a, b].real, T_out[key][i, a, b].imag)])
               for key in T_out for i, w in enumerate(om)])
    # the same T-matrix from the upfolded matrix (m = MMAX, eta = ETA) of the finite lattice
    pos = {int(c): n for n, c in enumerate(C)}
    sup = [pos[(L // 2 + dx) * L + (L // 2 + dy)] for dx, dy in SUPPORT]
    zt = (TEST_OMEGAS + 1j * ETA) ** 2
    G0c = cluster_green_modes(*modes_max["host"], len(C), zt)
    G0inv = np.linalg.inv(G0c)
    gt = host_green(TEST_OMEGAS, ETA)
    for key, _, Md, Kd in DEFECTS:
        mvec = np.ones(len(C))
        mvec[i0] = Md
        gu = cluster_green_modes(*modes_max[key], len(C), zt) / np.sqrt(np.outer(mvec, mvec))
        Tc = (G0inv @ (gu - G0c) @ G0inv)[:, sup][:, :, sup]
        Td, _ = defect_tmatrix(TEST_OMEGAS, ETA, Md, Kd, gt)
        dev = np.linalg.norm(Tc - Td, axis=(1, 2)) / np.linalg.norm(Td, axis=(1, 2))
        print(f"T-matrix {key}: upfolded (m={MMAX}, eta={ETA}) vs direct, max relative deviation {dev.max():.1e}")
    print(f"spectral figure: n_d = {ND}, eta = {ETA_A}, first order above w = {w_cut:.3f}; heavy resonance of one defect at {w_res:.3f}")
    print("wrote site/results/toy-local-dos.*, toy-convergence.*, toy-dispersion.*, toy-spectral.*, site/data/toy-*.csv")


if __name__ == "__main__":
    main()
