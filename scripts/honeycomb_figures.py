#!/usr/bin/env python3
"""Figures and data for the honeycomb toy page (run: npm run honeycomb, a few seconds).

Writes site/results/honeycomb-dispersion.{svg,pdf,png} (host dispersion along Γ-M-K-Γ beside the
density of states, both split into in-plane and out-of-plane motion) and the plotted numbers as
CSV under site/data/. Same style and colours as scripts/toy_figures.py.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from phdef import honeycomb as hc  # noqa: E402
from toy_figures import plt, INK, INK2, MUTED, AXIS, save, write_csv  # noqa: E402

COLOR = {"in": "#2a78d6", "out": "#eb6834"}          # categorical slots 1-2, validated
NGRID, SIGMA = 300, 0.006


def branches(k):
    """In-plane (4) and out-of-plane (2) frequencies at the k-points k (n, 2)."""
    D = np.array([hc.bloch(kk) for kk in k])
    w = lambda idx: np.sqrt(np.clip(np.linalg.eigvalsh(D[:, idx][:, :, idx]), 0.0, None))
    return w(hc.IN_PLANE), w(hc.OUT_OF_PLANE)


def main():
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


if __name__ == "__main__":
    main()
