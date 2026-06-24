# Electron–defect $T$-matrix: diagonal $T$, self-energy & spectral function (MoS₂ V$_S$, O$_S$)

Beyond-Born electron–defect scattering for monolayer MoS₂, computed by **explicit-summation downfolding** + **Wannier interpolation** — the production path of the Anvil EDT project, ported to Kestrel and validated against it to the meV. From the same EDI matrix element $M_{n\mathbf k,m\mathbf k'}=\langle n\mathbf k|\Delta V|m\mathbf k'\rangle$ used for the [defect-DOS diagonalization](vs-defect-dos-diag.html), we resum **all orders** of multiple scattering and obtain the on-shell diagonal $T$-matrix, the defect self-energy, and the spectral function $A(k,\omega)$ along $\Gamma$–M–K–$\Gamma$.

## 1. Method — explicit-summation downfolding + Wannier interpolation

The bare $M$ over $\sim$60 bands is split into an **active** block $P$ (the 11 Wannier bands 7–17 spanning the gap region) and a **rest** block $Q$ (bands 18–70). The rest is folded into an *exact* energy-dependent self-energy on the active space (Feshbach), evaluated statically at a reference energy $\omega_0=-5.955$ eV (near the VBM):

$$\Sigma(\omega_0)=W_{PQ}\,(\omega_0-E_Q-W_{QQ})^{-1}\,W_{QP},\qquad W=M\,\mathrm{Ry}/N_{\mathbf k}$$

computed by **one** eigendecomposition of $H_{QQ}=E_Q+W_{QQ}$ (no Krylov/Sternheimer — that route stagnates at the gap). This gives the dressed 11-band potential $\tilde V = M_{PP}+N_{\mathbf k}\Sigma(\omega_0)$. The active block is then Wannier-interpolated to arbitrary $\mathbf k$ and the multiple scattering resummed exactly by a small inversion (Koster–Slater):

$$T(\mathbf k,\mathbf k;\omega)=\big[\,1-\tilde V\,G^{A}(\omega)\,\big]^{-1}\tilde V,\qquad
\Sigma_{\rm defect}(\mathbf k,\omega)=n_d\,T(\mathbf k,\mathbf k;\omega),\qquad
A(\mathbf k,\omega)=-\frac{1}{\pi}\,\mathrm{Im}\,\mathrm{Tr}\,\big[\omega-H_0(\mathbf k)-n_d T\big]^{-1}.$$

**Validation (exact vs Anvil).** The downfolded in-gap levels reproduce the independent EDT code: V$_S$ $+1.209$ eV (Anvil explicit-60 $+1.205$, DFT $+1.19$); O$_S$ $+0.730/+0.725$ eV (Anvil unrelaxed $+0.731/+0.725$ — our O$_S$ is frozen geometry). The rest-space dressing moves the *bare* 11-band level (V$_S$ $+1.484$) down to the band-converged value — the beyond-Born content this captures.

## 2. Diagonal $T$-matrix along $\Gamma$–M–K–$\Gamma$

On-shell $T(\mathbf k,\mathbf k;\varepsilon_{\rm top}(\mathbf k))$ for the top valence band, **with** rest-space ($\tilde V$, solid) vs **without** (bare $M$, dashed).

![V_S diagonal T-matrix](../assets/tmat_tpath_vs.png)
*MoS₂ **V$_S$**: the $T$-matrix peaks sharply at **K** (the VBM), Im $T_{PP}\approx-0.065$ Ry. The rest-space curves (solid) are systematically stronger than bare (dashed) — V$_S$ carries a **strong** rest dressing.*

![O_S diagonal T-matrix](../assets/tmat_tpath_os.png)
*MoS₂ **O$_S$**: $\sim$3–4× stronger scattering (Im $T_{PP}\approx-0.24$ Ry at K), but solid$\approx$dashed — the rest dressing is **mild** for the isovalent substitution.*

## 3. Defect self-energy $n_d\Sigma_{\rm VBM}(\mathbf k)$

On-shell self-energy of the VBM band at defect density $n_d=2.78\%$ (one defect per $6\times6$ cell): **Re = level shift**, **Im = defect-limited broadening**.

![V_S self-energy](../assets/tmat_selfe_vs.png)
*V$_S$: VBM broadening peaks near K at $\sim$31 meV (shift $\sim$38 meV); elsewhere only a few meV.*

![O_S self-energy](../assets/tmat_selfe_os.png)
*O$_S$: much stronger — broadening up to $\sim$90 meV near K (shift $\sim$55 meV). The K-point VBM is most affected because that is where the defect states derive from.*

## 4. Spectral function $A(\mathbf k,\omega)$

$A(\mathbf k,\omega)=-\frac{1}{\pi}\mathrm{Im}\,\mathrm{Tr}\,[\omega-H_0-n_d T]^{-1}$ along $\Gamma$–M–K–$\Gamma$, **with** vs **without** rest-space (log color; cyan dashed = bare $\varepsilon_{n\mathbf k}$).

![V_S spectral function](../assets/tmat_spectral_vs.png)
*V$_S$: host bands acquire a $\mathbf k$-dependent defect linewidth, plus a flat **in-gap defect resonance** ($\sim$VBM$+1.2$ eV).*

![O_S spectral function](../assets/tmat_spectral_os.png)
*O$_S$: visibly heavier band-edge smearing near K, tracking its larger $T$-matrix.*

## 5. Summary

| Quantity | V$_S$ (vacancy) | O$_S$ (substitution) |
|---|---|---|
| Downfolded in-gap level (vs Anvil) | $+1.21$ eV ✓ ($+1.205$) | $+0.73$ eV ✓ ($+0.731$) |
| Im $T_{PP}$ near K | $\approx-0.065$ Ry | $\approx-0.24$ Ry (~3.7×) |
| VBM broadening near K ($n_d=2.78\%$) | $\sim$31 meV | $\sim$90 meV |
| Rest-space dressing | **strong** ($\sim$1.7×) | mild |

**Headline:** the isovalent **O$_S$ scatters the valence-band edge far harder than the V$_S$ vacancy**, while **V$_S$'s scattering depends much more on the rest-space dressing** — precisely the beyond-Born physics the explicit-summation $T$-matrix resolves. Pipeline and matrix elements are in place to repeat this for WS₂ V$_S$/O$_S$.
