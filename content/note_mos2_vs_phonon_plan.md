# Phonon spectral function of MoS₂ with S vacancies — computational plan

Goal: compute $B_{\mathbf q\nu}(\omega)$ for monolayer MoS₂ containing a dilute concentration of sulfur vacancies (V$_{\rm S}$), using the defect T-matrix formalism derived in the companion note ([FM self-energy & defect phonons](fm-phonon-defect.html), §4–§8): host resolvent $g(z)$ from primitive-cell DFPT, defect perturbation $V(z)=\Delta\mathcal D-z\,\varepsilon$ from a supercell force-constant difference, exact small-block T-matrix, dilute configurational average $\pi_{\mathbf q\nu}=c\,t_{\mathbf q\nu}$, spectral function from the boxed formula of §6. This page is the plan + TODO list; results will get their own page.

## 1. Ingredients already in hand

The electron-side EDI study of the *same defect* (MoS₂ + V$_{\rm S}$, validated against the official reference) fixed the structural model and the DFT settings; we reuse them verbatim for maximum cancellation in $\Delta\Phi$:

| Ingredient | Value / location | Status |
|---|---|---|
| Code | QE 7.5 at `~/qe-7.5` (pw.x, **ph.x**, q2r.x, matdyn.x built) | ready |
| Primitive cell | $a=3.18518$ Å, 24 Å vacuum, `assume_isolated='2D'` | ready (EDI run) |
| Pristine supercell | $6\times6$, 108 atoms, $\Gamma$-only | ready (EDI run) |
| Defect supercell | $6\times6$ + V$_{\rm S}$, 107 atoms, $\Gamma$-only, smearing 0.01 Ry | geometry to re-verify (relaxation!) |
| Pseudos / cutoff | PseudoDojo NC v0.5 PBE stringent, $E_{\rm cut}=100$ Ry | ready |
| Defect concentration | $n_d=10^{12}\,{\rm cm^{-2}}$ (matching the EDI transport run) $\Rightarrow c=n_d A_{\rm cell}\approx8.8\times10^{-4}$/cell | fixed |
| Electron–defect rates (for the optional e–ph bubble) | EDI `mos2_inv_tau.dat` | ready |

A $\Gamma$-only $6\times6$ supercell samples force constants exactly on a $6\times6$ commensurate grid — so the host DFPT $\mathbf q$-grid must be $6\times6\times1$ to make $\Delta\Phi$ well defined site-by-site.

## 2. Workflow

### P0 — Host phonons (primitive cell, DFPT)

ph.x on the relaxed primitive cell: $6\times6\times1$ $\mathbf q$-grid, with the 2D Coulomb treatment consistent with the SCF (`assume_isolated='2D'`; Born charges + 2D LO–TO via `loto_2d` in q2r.x/matdyn.x — monolayer MoS₂ is polar). Then q2r.x → real-space IFCs $\Phi^0$; matdyn.x along $\Gamma$–M–K–$\Gamma$ and on dense grids. *Output*: $\Phi^0$, $\omega_{\mathbf q\nu}$, $\mathbf e_{\mathbf q\nu}$, Fourier-interpolable everywhere.

### P1 — Defect supercell: relax, then finite-displacement force constants

1. **Relax** the 107-atom defect supercell (BFGS, same cutoff/pseudos/smearing; forces < 10⁻⁴ Ry/Bohr). Phonons are second derivatives at the *minimum* — the EDI run only needed the SCF potential, so the shipped geometry's force residual must be checked first.
2. **Finite displacements** (phonopy with the QE-pw.x interface, displacement 0.015 Å): one forces-only SCF per irreducible displacement. $C_{3v}$ symmetry of V$_{\rm S}$ reduces the 321 displacements to an irreducible set (expect order 10²). Same for the **pristine** 108-atom supercell (full crystal symmetry → only a handful of displacements); computing both supercells with the *same method and settings* makes systematic errors cancel in the difference.

*Output*: $\Phi^{\rm def}$ (107 atoms), $\Phi^{\rm prist}$ (108 atoms).

### P2 — Assemble the defect perturbation $V(z)$

Map defect-supercell sites onto pristine sites (vacancy site has no image); build $\Delta\Phi=\Phi^{\rm def}-\Phi^{\rm prist}$ on the common 106 atoms + the vacancy row/column handled by projection (the vacant S contributes $\Delta M=-M_{\rm S}$ and decoupled rows — project its 3 degrees of freedom out of the cluster space). Then: (i) plot $\lVert\Delta\Phi(d)\rVert$ vs distance $d$ from the vacancy — require ≥ 2 orders of decay before the supercell boundary; (ii) truncate to a defect cluster (vacancy + Mo first shell + S shells out to $\sim5$–6 Å, ~25–40 atoms ⇒ cluster dimension ~75–120); (iii) re-impose the acoustic sum rule on the truncated $\Delta\Phi$. Build $\varepsilon$ (mass term) and $V(z)=\Delta\mathcal D-z\,\varepsilon$.

### P3 — Host resolvent on the cluster (spectral density + Hilbert transform)

The T-matrix inversion lives in the real-space cluster basis (that is what makes it a small block — see P4), so only the $d\times d$ cluster block $g_{ab}(z)$ is ever needed. Here $a\equiv(l,\kappa,\alpha)$ is a composite site index — the (host-mass-scaled) displacement of atom $\kappa$ in cell $l$ along Cartesian direction $\alpha$, restricted to the truncated defect cluster (vacancy + neighbour shells, $d=3N_{\rm cluster}\approx75$–$120$ after projecting out the 3 vacant-site degrees of freedom). It is the same orthonormal basis in which $V(z)=\Delta\mathcal D-z\,\varepsilon$ is expressed, with mode projections $\langle l\kappa\alpha|\mathbf q\nu\rangle=e_{\nu\alpha}(\kappa;\mathbf q)\,e^{i\mathbf q\cdot\mathbf R_l}/\sqrt{N_p}$ (cf. §4 of the [derivation note](fm-phonon-defect.html)). Compute it in **two steps**, not as one BZ sum per frequency:

1. **One pass** over a dense Fourier-interpolated $\mathbf q$-grid (start $300\times300\times1$, converge; tetrahedron or adaptive smearing) accumulates the cluster-projected spectral density on a fine $\lambda=\omega'^2$ mesh:

$$
\rho_{ab}(\lambda)=\frac{1}{N_p}\sum_{\mathbf q\nu}\langle a|\mathbf q\nu\rangle\langle\mathbf q\nu|b\rangle\,N_p\;
\delta\big(\lambda-\omega_{\mathbf q\nu}^2\big).
$$

2. The resolvent at **any** complex $z$ then follows by a Hilbert (Kramers–Kronig) transform,

$$
g_{ab}(z)=\int_0^{\lambda_{\max}}\! d\lambda\;\frac{\rho_{ab}(\lambda)}{z-\lambda},
\qquad
\mathrm{Im}\,g_{ab}\big((\omega+i0^+)^2\big)=-\pi\,\rho_{ab}(\omega^2),
$$

so the scan over thousands of $z$ points for the T-matrix is essentially free, and the small-$\eta$ convergence burden moves onto the (cheap) $\lambda$-mesh of $\rho$ instead of the $\mathbf q$-sum. Built-in sanity checks: completeness $\int d\lambda\,\rho_{ab}(\lambda)=\delta_{ab}$, positivity $\rho_{aa}(\lambda)\ge0$, and Hermiticity of $g$. Frequency mesh for the downstream T-matrix: 0–65 meV (top of optical branches ≈ 58 meV), step ≤ 0.05 meV near resonances.

### P4 — T-matrix and self-energy

Per frequency: invert $[\mathbb 1-V(z)g(z)]$ on the cluster, $T(z)=[\mathbb 1-Vg]^{-1}V$; project onto host modes $\to t_{\mathbf q\nu}(z)$; $\pi_{\mathbf q\nu}=c\,t_{\mathbf q\nu}$ with $c=8.8\times10^{-4}$ (and a concentration series $10^{11}$–$10^{13}\,{\rm cm^{-2}}$ for trends).

### P5 — Observables

$B_{\mathbf q\nu}(\omega)$ maps along $\Gamma$–M–K–$\Gamma$ (color maps, host dispersion overlaid); mode linewidths $\Gamma_{\mathbf q\nu}$ and lifetimes $1/\tau=2\Gamma$ vs $\omega$ (input for defect-limited lattice thermal conductivity); resonant-mode search $\det[\mathbb 1-\mathrm{Re}\,V g]=0$ (expect vacancy-induced quasi-localized modes in/below the acoustic–optical gap); $\Delta$DOS from the Krein determinant.

### P6 — Optional (second pass): non-adiabatic e–ph bubble

Add $2\omega_{\mathbf q\nu}[\Pi^{ep}(\omega)-\Pi^{ep}(0)]$ with electron lines broadened by the EDI electron–defect rates. Needs $g_{mn\nu}(\mathbf k,\mathbf q)$ (EPW is built) — deferred until the structural part is validated; at $n_d=10^{12}\,{\rm cm^{-2}}$ the carrier density is low and this term should be a small correction.

## 3. Validation gates

| Gate | Check | Pass criterion |
|---|---|---|
| V0 | Pristine primitive dispersion (P0) | no imaginary modes; quadratic ZA near $\Gamma$ (rotational ASR); branches match published monolayer-MoS₂ phonons |
| V1 | Pristine supercell $\Gamma$ phonons (phonopy, P1) vs folded DFPT (P0) | mode-by-mode agreement to a few cm⁻¹ — validates the finite-displacement machinery before trusting $\Phi^{\rm def}$ |
| V2 | $\Delta\Phi$ decay + ASR residual (P2) | ≥ 2 orders down at the cluster edge; ASR residual < 10⁻⁴ of NN force constant |
| V3 | Limits of the T-matrix code (P4) | $V\to0$ or $c\to0$: $B\to\delta(\omega-\omega_{\mathbf q\nu})$, $\int B\,d\omega=1$ per mode |
| V4 | Born + mass-only artificial test (³⁴S isotope on the S site) | numerical $\Gamma_{\mathbf q\nu}$ reproduces the analytic Tamura formula |
| V5 | Defect supercell direct diagonalization (321 modes of $\Phi^{\rm def}$) | resonant/localized mode energies from the T-matrix match the supercell fingerprint; Krein mode count conserved |

## 4. Cost estimate (Kestrel, 1 shared node = 104 cores; supercell SCF ≈ 5 min from the EDI run)

| Step | Work | Estimate |
|---|---|---|
| P1.1 defect relax | ~20–40 BFGS steps | 2–4 node-h |
| P0 DFPT primitive | ~7 irreducible q × 9 modes | 2–6 node-h |
| P1.2 phonopy pristine | ~4–10 displacements | < 1 node-h |
| P1.2 phonopy defect | ~100–200 displacements × 5 min | 8–17 node-h (job-array friendly) |
| P2–P5 post-processing | Python/numpy on login node | cheap |

Total ≈ 15–30 node-hours on `shared` — comparable to the EDI electron-side study.

## 5. TODO

- [x] **T0** Install phonopy (`pip install --user phonopy`) and check the pw.x interface on the primitive cell
- [x] **T1** (P1.1) Check force residual of the shipped 107-atom geometry; relax the defect supercell
- [x] **T2** (P0) DFPT: ph.x $6\times6\times1$ on the primitive cell with 2D polar treatment; q2r + matdyn; gate **V0**
- [x] **T3** (P1.2) phonopy finite displacements: pristine 108-atom supercell; gate **V1** vs folded DFPT
- [x] **T4** (P1.2) phonopy finite displacements: relaxed defect supercell (job array)
- [x] **T5** (P2) `tools/build_dphi.py` — site mapping, $\Delta\Phi$, decay plot, cluster truncation, ASR; gate **V2**
- [x] **T6** (P3) `tools/host_resolvent.py` — single-pass cluster spectral density $\rho_{ab}(\lambda)$ (tetrahedron) + Hilbert transform to $g_{ab}(z)$; checks: completeness $\int\rho_{ab}\,d\lambda=\delta_{ab}$, positivity, $\lambda$-mesh / $\mathbf q$-grid convergence
- [x] **T7** (P4) `tools/tmatrix_phonon.py` — $V(z)$, block inversion, mode projection, $\pi_{\mathbf q\nu}$; gates **V3, V4**
- [x] **T8** (P5) `tools/spectral.py` — $B(\mathbf q,\omega)$ maps, linewidths, resonant modes, $\Delta$DOS; gate **V5**
- [x] **T9** Results page: spectral-function maps + linewidths + resonance table, uploaded here
- [ ] **T10** *(optional, P6)* Non-adiabatic e–ph bubble with EDI-dressed electron lines

## 6. Risks and open questions

1. **In-gap electronic states of V$_{\rm S}$ + smearing.** The vacancy puts deep states in the gap; with 0.01 Ry smearing their fractional occupation can contaminate forces. Check the occupation spectrum after relaxation; if states sit at $E_F$, compare forces against `occupations='fixed'` (the neutral cell should be insulating).
2. **2D polar corrections.** The host IFC interpolation must use the 2D (not 3D) dipole–dipole treatment, or the LO branch near $\Gamma$ — and hence $g(z)$ — is wrong. The *difference* $\Delta\Phi$ of the neutral vacancy should be short-ranged, but gate V2 is the explicit check.
3. **Supercell size.** If V2 fails (slow $\Delta\Phi$ decay), a $8\times8$ defect supercell is the fallback — ~3× cost on P1.
4. **ZA branch numerics.** The quadratic flexural branch is fragile: enforce the rotational-invariance (Huang) sum rules in matdyn/post-processing, and watch $g(z)$ convergence at small $\omega$ where the ZA DOS diverges as a 2D van Hove tail.
