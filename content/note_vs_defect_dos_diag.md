# Defect DOS by single-defect Hamiltonian diagonalization (MoS₂ V$_S$)

The literature electron–defect approach (e.g. eq. (S1) of the referenced works) writes a single-defect Hamiltonian
$$
\hat H_{e\text{-}i}=\sum_{ij}c_i^\dagger c_j\,g_{ij},\qquad g_{ij}=\langle i|\Delta V|j\rangle ,
$$
and obtains the defect density of states by diagonalizing it in the coarse-grid Bloch basis $i=(n,\mathbf k)$:
$$
H_{(n\mathbf k),(m\mathbf k')}=\varepsilon_{n\mathbf k}\,\delta+g_{(n\mathbf k),(m\mathbf k')},\qquad
g_{ij}=M_{ij}\,\frac{\mathrm{Ry}}{N_{\mathbf k}} .
$$
On an $N_{\mathbf k}$-point coarse grid this is **equivalent to a single defect in an $N_{\mathbf k}$-cell supercell**, obtained from one dense diagonalization (Koster–Slater / explicit $T$-matrix). $\Delta V$ and $M=\langle n\mathbf k|\Delta V|m\mathbf k'\rangle$ come from the EDI code (direct mode); $\varepsilon_{n\mathbf k}$ are the pristine host bands.

> **Two corrections got this page to its current (validated) state.** Both were found by cross-checking against an independent implementation (the EDT/Sternheimer code on a separate machine).
> 1. **Normalization** $g=M\,\mathrm{Ry}/N_{\mathbf k}$ — the stored $M$ is per-cell and in Rydberg; the single-defect coupling carries the $1/N_{\mathbf k}$ defect-density factor and the Ry→eV conversion. Omitting it over-weighted the perturbation ${\sim}10.6\times$.
> 2. **A bug in EDI direct mode** left the *bra* wavefunction stale (frozen at the last $k$-point) for every $\langle n\mathbf k|\Delta V|m\mathbf k'\rangle$, making the matrix $M$ **non-Hermitian** ($\lVert M-M^\dagger\rVert/\lVert M\rVert\approx2$) — which broke the crystal $C_3$ symmetry and **split the $e$-doublet**. Fixed (see §1).
>
> All earlier numerical tables on this page (from before these fixes) are retracted; the result below is the corrected, cross-validated one.

## 1. The Hermiticity bug and its fix

For the same $\mathbf k$ block, $\Delta V$ is Hermitian so $M(n,m)$ **must** equal $M(m,n)^\ast$. It did not. The breaking was localized step by step:

| object | $C_3$-symmetric? |
|---|---|
| $\Delta V$ cube ($V_d-V_p$, supercell grid) | ✅ exact (QE symmetrizes $V$) |
| folded $V_{\rm folded}$ (supercell → primitive) | ✅ exact |
| host bands $\varepsilon_{n\mathbf k}$ | ✅ |
| **EDI-direct $M$** | ❌ broken (both local and nonlocal parts) |

The same-$\mathbf k$ block $M(\mathbf k)$ was non-Hermitian (rel $\approx2.0$, near-anti-Hermitian). Root cause: in `ed_direct_from_files`, the bra wavefunctions (`psir_ki`/`becd_ki`/`becpc_ki`) were broadcast in a pre-loop that overwrote a single buffer, leaving only the **last** $k$-point; the compute loop then used that stale bra for **all** $(\mathbf k_i,\mathbf k_f)$. So every $M(\mathbf k_i,\mathbf k_f)$ used bra $=\psi_{\mathbf k_{N}}$. (The independent EDT code reads the bra fresh for each $\mathbf k_i$ — correct.) **Fix:** process only pairs whose $\mathbf k_i$ is local to the MPI pool, taking the bra straight from the cache (no broadcast) and the ket from the streamed panel — zero extra communication, and the edmat is actually faster.

## 2. Corrected V$_S$ defect DOS

With the bra-fixed, Hermitian, $C_3$-symmetric $M$ (full 1–66 band window, 12×12 grid = 144-cell):

![V_S defect DOS (bra-fixed)](../assets/dos_vs_brafixed.png)
*Host vs defect DOS (left) and Friedel ΔDOS (right), referenced to the VBM. In-gap defect states: an $a_1$ singlet and a **degenerate $e$-doublet** (the taller in-gap peak = 2 states).*

| in-gap state | $E-E_{\rm VBM}$ (eV) |
|---|---|
| $a_1$ (singlet) | +0.24 |
| $e$ (doublet, **degenerate**) | **+1.196, +1.196** |

The $e$-doublet is degenerate (split $=0.000$ eV) and lands at +1.20 eV — in the DFT-supercell range (9×9 ≈ 1.06; 6×6 ≈ 1.15) and matching the independent EDT value below.

## 3. Cross-validation against an independent code (EDT)

The same coarse Bloch $M$ computed by the independent EDT implementation, and the bra-fixed EDI $M$, agree on the symmetry that matters:

| check | EDI before fix | EDI after fix | EDT (independent) |
|---|---|---|---|
| same-$\mathbf k$ block Hermiticity $\lVert M-M^\dagger\rVert/\lVert M\rVert$ | 2.0 | $\sim10^{-13}$ | $\sim10^{-14}$ |
| $M$-block spectra form $C_3$ triplets (of 144 $k$) | 0/144 | **143/144** | 143/144 |
| $e$-doublet | split 0.4 eV | **+1.196 ×2** | +1.205 ×2 |

The $e$-doublet now agrees with EDT to ~0.01 eV. The residual $a_1$ offset (+0.24 here vs +0.005 in EDT) is the alignment choice — these runs use `pot_align='none'`, EDT uses `'vacuum'` (a rigid constant shift; it does not affect the degeneracy or the $e$-level).

## 4. Setup & caveats

```
edmat_direct_from_file = .true.   ! direct-mode Bloch M, band_ed='1-66'
filki_direct = filkf_direct = 'kfull.dat'   ! 144 coarse k (12×12)
pot_align = 'none'
```
Offline: $H=\mathrm{diag}(\varepsilon)+M\,\mathrm{Ry}/N_{\mathbf k}$, `eigvalsh`; in-gap = eigenvalues between host VBM and CBM.

- **Frozen, unrelaxed** neutral-$V_S$ $\Delta V$; absolute levels carry the usual PBE / Γ-only-SCF / coarse-$k$ uncertainties.
- `pot_align='none'` shifts the $a_1$ vs an aligned reference (does not change the $e$-doublet).
- The O$_S$ and Si-vacancy controls are being re-run with the bra-fixed code; results to follow.
