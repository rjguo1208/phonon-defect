# Defect DOS by single-defect Hamiltonian diagonalization (MoS₂ V$_S$)

The literature electron–defect approach (e.g. eq. (S1) of the referenced works) writes a single-defect Hamiltonian
$$
\hat H_{e\text{-}i}=\sum_{ij}c_i^\dagger c_j\,g_{ij},\qquad g_{ij}=\langle i|\Delta V|j\rangle ,
$$
and obtains the defect density of states by diagonalizing it. Here we realize this for the S vacancy in monolayer MoS₂, reusing the difference potential $\Delta V$ and the Bloch matrix elements $g_{ij}$ already produced by the EDI code. The object actually diagonalized is the **total** Hamiltonian in the coarse-grid Bloch basis $i=(n,\mathbf k)$,
$$
H_{(n\mathbf k),(m\mathbf k')}=\varepsilon_{n\mathbf k}\,\delta+g_{(n\mathbf k),(m\mathbf k')},\qquad
\rho(E)=-\tfrac1\pi\,\mathrm{Im\,Tr}\,(E+i\eta-H)^{-1},
$$
with $g$ taken from EDI's `mos2_edmat_bloch.dat` ("Bloch M before Wannier rotation") and $\varepsilon_{n\mathbf k}$ from the primitive nscf. On the coarse $12\times12$ k-grid this is **equivalent to a single defect in a 144-cell ($12\times12$) supercell** — larger than the explicit $6\times6$/$9\times9$ DFT supercells, but obtained at the cost of one $1584\times1584$ diagonalization (seconds), because $\Delta V$ and the host bands are reused rather than recomputed self-consistently.

## 1. Method ⇄ our other work

This is the DOS form of the **Koster–Slater / explicit T-matrix** approach: $g_{ij}=\langle n\mathbf k|\Delta V|m\mathbf k'\rangle$ is exactly the explicit Bloch $M$-matrix used on the [electron side](https://rjguo1208.github.io/) (claude-sternheimer), and the electronic mirror of the phonon [defect T-matrix](fm-phonon-defect.html). The single defect breaks translational symmetry, so $g$ is dense in $(\mathbf k,\mathbf k')$; diagonalizing it folds back to the $12\times12$ supercell levels.

## 2. Band-space convergence: 5-band vs 11-band Wannier

The defect $e$ states are notoriously sensitive to the size of the band/Wannier space — too small a manifold **over-screens** them (drags them toward the VBM). We compare two Wannier setups built from the same $\Delta V$:

| | 5-band (Mo:d, bands 13–17) | 11-band (Mo:d + S:p, bands 7–17) |
|---|---|---|
| in-gap defect levels ($E-E_{\rm VBM}$, eV) | 0.06, 0.09, **0.50, 0.50, 0.60, 0.60** | 0.08, 0.19, 0.19, **0.82, 0.82, 0.83** |
| character | a₁ near VBM + over-screened $e$ | a₁ + upward-shifted $e$ manifold |

![5-band defect DOS](../assets/defect_dos_5band.png)
*5-band (Mo:d): host vs defect DOS (left), Friedel ΔDOS (right). Purple = in-gap defect levels. The $e$ manifold sits at ~0.5–0.6 eV (over-screened).*

![11-band defect DOS](../assets/defect_dos_11band.png)
*11-band (Mo:d + S:p): the in-gap $e$ manifold has moved up to ~0.82 eV.*

**Trend.** Enlarging the manifold from 5 to 11 bands (adding S:p and the lower valence bands) pushes the $e$ states from ~0.5–0.6 eV up to ~0.82 eV — **partially correcting the over-screening**, and moving toward the explicit/DFT references:

| reference | $e$ level ($E-E_{\rm VBM}$) |
|---|---|
| 5-band diag | 0.50–0.60 eV |
| **11-band diag (this work)** | **0.82 eV** |
| 9×9 DFT supercell | ~1.06 eV |
| 6×6 DFT supercell | ~1.15 eV |
| explicit 21-band $T$-matrix (claude-sternheimer) | ~1.35 eV |

The monotonic climb with band-space size is consistent with the claude-sternheimer finding that the deep $e$ doublet requires a large explicit manifold (≈21 bands) to fully escape over-screening. 11 bands gets roughly halfway from the 5-band result to the DFT value — a clear, controlled improvement, not yet full convergence.

## 3. Practical notes (why this was expensive to set up)

Generating the 11-band $g_{ij}$ with EDI's `edi.x` hit four independent memory/MPI walls on Kestrel, each masking the next:

1. **Full-core OOM** — 104 ranks/node is too many; run 26 ranks/node (reserve the full node with `--ntasks-per-node=104 --mem=0`, launch fewer with `srun --ntasks-per-node=26`).
2. **Cray MPICH CMA crash** (`process_vm_readv copy size mismatch` on a ~253 MB transfer) → `export MPICH_SMP_SINGLE_COPY_MODE=NONE`.
3. **Fine-interpolation array** — `fine_nk=300` preallocates ~174 GB on one rank; the diagonalization needs only the coarse $M$, so `fine_nk=12`.
4. **psir cache** — wavefunctions are cached on the *supercell* $240^3$ cube grid (to match $\Delta V$): total $144\,k\times11\,{\rm band}\times240^3\times16\,{\rm B}=437$ GB, spread over nodes. A Kestrel standard node has 243 GB, so 2 nodes (219 GB/node) + overhead OOMs; **3 nodes (146 GB/node)** fit.

Winning configuration: 3 nodes, 26 ranks/node, `fine_nk=12`, single-copy off, full-node memory. The diagonalization itself is a trivial $1584\times1584$ eigenproblem afterwards.

## 4. Caveats

- $\Delta V$ is the frozen $6\times6$ supercell difference potential (neutral V$_S$, short-ranged — fine; charged defects would need range separation).
- The $g_{ij}$ here is the **bare** Bloch matrix element (before Wannier rotation); the Wannier step is used only to validate band interpolation, not for the diagonalization.
- Band-space (11) is still not fully converged for the deep $e$ states — the absolute level (~0.82 eV) is a lower bound relative to the DFT/explicit value; the *method* and the *over-screening-correction trend* are the robust results.
- Γ-only supercell SCF + smearing for $\Delta V$; coarse $12\times12$ Bloch grid for the diagonalization.

Artifacts: `/scratch/rjguo/edi_example/edi_run/edi/` — `mos2_edmat_bloch.dat` (11-band $g_{ij}$), `diagonalize_dos.py`, `defect_dos_{5,11}band.{png,npz}`.
