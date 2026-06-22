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

## 4. Band-manifold convergence (corrected)

Slicing the host-band basis and re-diagonalizing shows the $e$-doublet converging with manifold size — **with the bra-fixed $M$ the level descends toward the DFT value from above**:

| manifold | $e$-doublet ($E-E_{\rm VBM}$, eV) |
|---|---|
| bands 7–17 (11) | +1.49 |
| bands 7–66 (60) | +1.21 |
| bands 1–66 (66) | +1.196 |
| DFT supercell (9×9 / 6×6) | 1.06 / 1.15 |

A larger manifold gives more complete screening, **lowering** the level toward the DFT range. (An earlier version of this page reported the *opposite* trend — the level "climbing" 0.5→0.8→1.3 eV — which was an artifact of the non-Hermitian $M$; that is retracted.)

## 5. Defect-induced electron self-energy (T-matrix)

The same bra-fixed coupling $g$ gives the defect self-energy in the dilute (single-site) limit, $\Sigma(\omega)=n_d\,T(\omega)$, $T(\omega)=g\,[1-G_0(\omega)g]^{-1}$, $n_d=1/N_{\mathbf k}$ — the complex, frequency-dependent complement of the (real-spectrum) diagonalization.

![V_S T-matrix self-energy](../assets/selfenergy_vs_166.png)
*Tr $\Sigma(\omega)$ (left): dispersive Re $\Sigma$ and Im $\Sigma$ dips at **both** bound states. Friedel ΔDOS (right): two in-gap poles — the $a_1$ (+0.24) and the $e$-doublet (+1.2, ~2× higher = 2 states).*

The bound-state poles of $T$ ($\det[1-G_0g]=0$) coincide with the diagonalization eigenvalues, closing the two views — **both** the $a_1$ and the $e$-doublet:

| in-gap pole | diagonalization (DOS) | T-matrix (Friedel ΔDOS) | EDT |
|---|---|---|---|
| $a_1$ (singlet) | +0.238 | +0.243 | +0.005$^\dagger$ |
| $e$ (doublet) | +1.196 | +1.206 | +1.205 |

$^\dagger$the $a_1$ position is alignment-sensitive (`pot_align='none'` here vs `'vacuum'` in EDT — a rigid shift); the $e$-doublet is alignment-robust. So the **defect self-energy is obtained from the same $g$**: diagonalization gives the poles/DOS, the T-matrix gives the full complex $\Sigma(\omega)$ (level shift from Re $\Sigma$, scattering rate from Im $\Sigma$), and the two agree on both poles. (In a smaller band window, e.g. 7–66, the shallow $a_1$ is pushed to the VBM edge; the full 1–66 window resolves it cleanly at +0.24.)

## 6. Controls: O$_S$ (isovalent) and the Si vacancy

| system | DFT ground truth | single-defect diag (bra-fixed) | verdict |
|---|---|---|---|
| **V$_S$** | $e$ at 1.06–1.15 eV ([9×9 bands](supercell-bands.html)) | $a_1$ + $e$ **+1.196** | ✅ real defect states reproduced |
| **O$_S$** | **no in-gap state** ([9×9 bands](supercell-bands.html)) | spurious **+0.73** doublet | ❌ method artifact |
| Si vacancy (3D) | $t_2$-derived in-gap (supercell) | cluster near CBM (+0.6…+0.7) | ✅ in-gap manifold present |

![O_S defect DOS (bra-fixed)](../assets/dos_os_brafixed.png)
*O$_S$ defect DOS: only a weak spurious +0.73 in-gap feature (cf. the strong V$_S$ $e$-peak), while DFT has none.*

![Si vacancy defect DOS (bra-fixed)](../assets/dos_si_brafixed.png)
*Si vacancy (3D control): a $t_2$-derived in-gap manifold near the CBM — a genuine defect level present in the DFT supercell too.*

**The bra-fix makes V$_S$ correct, but O$_S$ still shows a spurious +0.73 in-gap doublet** while the dilute-limit DFT (9×9, even unrelaxed) has **none**. This is not a code bug, not the potential reference, and — as the band-character test below shows — **not the band-basis size either**. It is **defect-type-dependent**: whether the perturbation actually binds a state.

### Why O$_S$ fails: a real bound state vs a spurious band-edge hybrid

The diagonalization builds an in-gap level by letting $\Delta V$ **mix the pristine valence-band maximum (VBM) and conduction-band minimum (CBM)** into a state inside the gap. Decomposing the in-gap eigenvectors by band shows that **both** the V$_S$ $e$ and the O$_S$ "+0.73" peak at the band-edge states (band 13–14) with essentially **zero weight above band 30** — they are the *same kind* of low-band VBM–CBM hybrid. So the difference is **not** that O$_S$ "needs more bands" (an earlier, incorrect reading of mine): the 66-band basis represents both equally well.

The real distinction is whether that hybrid corresponds to a state that physically exists:

- **V$_S$ (vacancy)** removes an atom and leaves real dangling bonds; the VBM–CBM mixing produces a **genuine** bound state ($a_1$+$e$), confirmed by the DFT supercell. The diag reproduces it.
- **O$_S$ (isovalent substitution)** introduces no deep level: the *self-consistent* supercell redistributes charge and keeps the gap **clean**. But the diag applies the frozen $\Delta V$ as a **one-shot, static** perturbation to the *unperturbed* band edges and **over-mixes** VBM↔CBM, manufacturing an in-gap state the self-consistent DFT does not have.

So the static one-step treatment captures a real level when one exists (V$_S$, Si vacancy) but over-binds the band-edge coupling when none does (O$_S$). The checks below rule out every other candidate, leaving this defect-type mechanism.

- **Faithful elsewhere.** The 6×6 diag reproduces the DFT O$_S$ supercell band edges (VBM to ~10 meV, CBM to ~30 meV) and the valence bands to a few meV — it is faithful *everywhere except* the one spurious in-gap state.
- **Reference ruled out.** Primitive potential tiled to 6×6 vs the supercell pristine: RMS **1.8 meV** (~300× smaller than the O$_S$ $\Delta V$, RMS 0.58 eV).
- **Null/baseline test.** Pristine-vs-tiled-primitive ($\Delta V\approx0$) reproduces the host to **1.4 meV with zero in-gap states** — the method does not manufacture levels from nothing.
- **Sign test.** $-\Delta V$ does not remove the level, it relocates it (+0.73 → +1.0–1.15): a strong band-edge perturbation of *either* sign over-mixes.
- **$k$-grid test.** 12×12 → 6×6 leaves V$_S$ $e$ (+1.196→+1.191) and the O$_S$ doublet (+0.729→+0.729) fixed (figure below) — the level is set by the per-cell potential, not the supercell size.

![Convergence of in-gap levels to 90 bands: V_S vs O_S](../assets/conv90_viz.png)
*1–90 band convergence: both the V$_S$ $e$ and the O$_S$ doublet descend and stabilize. The O$_S$ doublet plateaus near +0.73 — it does not migrate out of the gap with more bands, consistent with a low-band band-edge hybrid (not a high-band effect).*

![Defect DOS vs k-grid: 6×6 vs 12×12](../assets/dos6_vs12.png)
*Defect DOS at 12×12 (solid) vs 6×6 (dashed): the in-gap peaks sit at the **same** energies for both grids — $k$-grid-independent.*

**Takeaway: the diagonalization cannot, on its own, tell a real defect level from a spurious band-edge hybrid.** Band convergence, $k$-grid, sign, and a clean baseline all look fine for V$_S$ (real) and O$_S$ (spurious) alike; only the DFT supercell ([9×9 bands](supercell-bands.html), O$_S$ = 0 in-gap) distinguishes them. The method is reliable for defects that bind a real level (V$_S$, Si vacancy) and over-binds isovalent perturbations that do not (O$_S$).

## 7. Setup & caveats

```
edmat_direct_from_file = .true.   ! direct-mode Bloch M, band_ed='1-66'
filki_direct = filkf_direct = 'kfull.dat'   ! 144 coarse k (12×12)
pot_align = 'none'
```
Offline: $H=\mathrm{diag}(\varepsilon)+M\,\mathrm{Ry}/N_{\mathbf k}$, `eigvalsh`; in-gap = eigenvalues between host VBM and CBM.

- **Frozen, unrelaxed** neutral-$V_S$ $\Delta V$; absolute levels carry the usual PBE / Γ-only-SCF / coarse-$k$ uncertainties.
- `pot_align='none'` shifts the $a_1$ vs an aligned reference (does not change the $e$-doublet).
- The O$_S$ and Si-vacancy controls are being re-run with the bra-fixed code; results to follow.
