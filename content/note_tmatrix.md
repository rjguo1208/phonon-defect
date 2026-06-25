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

**Gauge-locked Wannier rotation.** The interpolation requires the Wannier rotation $U(\mathbf k)$ to be in the *same gauge* as the Bloch states of the $M$ matrix. An initial run used a stale `filukk` (from a separate Wannierization) whose per-$\mathbf k$ phases did not match the 90-band NSCF — the real-space $\tilde V^W$ then failed to decay (Koster–Slater truncation captured only 55 %). Re-Wannierizing on the **same 90-band NSCF** (keeping exactly bands 7–17 → 11 WF, Mo:d + S:p) restores a gauge-consistent $U$; $\tilde V^W$ now decays $\sim$2 orders over $\sim$1.5 Å and $R_{\rm cut}{=}4$ captures $\geq$99 % (§6). All results below use the gauge-locked rotation.

## 2. Diagonal $T$-matrix along $\Gamma$–M–K–$\Gamma$

On-shell $T(\mathbf k,\mathbf k;\varepsilon_{\rm top}(\mathbf k))$ for the top valence band, **with** rest-space ($\tilde V$, solid) vs **without** (bare $M$, dashed).

![V_S diagonal T-matrix](../assets/tmat_tpath_vs.png)
*MoS₂ **V$_S$**: smooth, path-symmetric $T(k)$ (a hallmark of the gauge-locked, localized $\tilde V^W$). $|T_{PP}|$ reaches $\approx0.035$ Ry; at **K** Re $T_{PP}\to0$. With-rest (solid) and bare (dashed) differ in sign of Re at K ($-0.003$ vs $+0.011$ Ry) — V$_S$ carries a **non-trivial** rest dressing.*

![O_S diagonal T-matrix](../assets/tmat_tpath_os.png)
*MoS₂ **O$_S$**: stronger scattering, $|T_{PP}|$ up to $\approx0.080$ Ry and peaking near **K**; at the K-VBM $|T|\approx0.080$ vs V$_S$ $\approx0.004$ (**$\sim$19×**). Solid$\approx$dashed — the rest dressing is **mild** for the isovalent substitution.*

## 3. Defect self-energy $n_d\Sigma_{\rm VBM}(\mathbf k)$

On-shell self-energy of the VBM band at defect density $n_d=2.78\%$ (one defect per $6\times6$ cell): **Re = level shift**, **Im = defect-limited broadening**.

![V_S self-energy](../assets/tmat_selfe_vs.png)
*V$_S$: VBM broadening $|{\rm Im}\,n_d\Sigma|$ up to $\sim$9.5 meV; the level shift Re $n_d\Sigma$ ranges $-1$ to $+13$ meV across the path.*

![O_S self-energy](../assets/tmat_selfe_os.png)
*O$_S$: broadening up to $\sim$13 meV, and a markedly larger, uniformly negative **level shift** Re $n_d\Sigma\approx-16$ to $-30$ meV — the isovalent O pushes the VBM down harder than it broadens it.*

## 4. Spectral function $A(\mathbf k,\omega)$

$A(\mathbf k,\omega)=-\frac{1}{\pi}\mathrm{Im}\,\mathrm{Tr}\,[\omega-H_0-n_d T]^{-1}$ along $\Gamma$–M–K–$\Gamma$, **with** vs **without** rest-space (log color; cyan dashed = bare $\varepsilon_{n\mathbf k}$).

![V_S spectral function](../assets/tmat_spectral_vs.png)
*V$_S$: host bands acquire a $\mathbf k$-dependent defect linewidth, plus a flat **in-gap defect resonance** ($\sim$VBM$+1.2$ eV).*

![O_S spectral function](../assets/tmat_spectral_os.png)
*O$_S$: visibly heavier band-edge smearing near K, tracking its larger $T$-matrix.*

## 5. $T(nk,\omega)$ spectral map

The full energy dependence of the VBM-band diagonal $T(nk,\omega)$ along $\Gamma$–M–K–$\Gamma$ — **Re** (level shift) and **Im** ($-$Im $\propto$ scattering rate) as separate maps; dashed = on-shell $\varepsilon_{\rm VBM}(k)$.

![V_S T(k,omega) map](../assets/tmat_map_vs.png)
*V$_S$: a **dispersionless (flat-in-$k$) in-gap resonance** at $\omega\approx-4.9$ eV — Re $T$ shows the level-shift sign flip across it, Im $T$ a bright scattering ridge. Flatness in $k$ is the signature of a localized defect state (and a direct check that the gauge-locked interpolation is clean).*

![O_S T(k,omega) map](../assets/tmat_map_os.png)
*O$_S$: in-gap structure that is stronger near **K** and carries an additional branch toward the CBM, consistent with the larger O$_S$ $T$-matrix.*

## 6. Real-space locality & gauge validation

The Wannier interpolation and the $R_{\rm cut}$ truncation are only valid if the dressed potential $\tilde V^W(\mathbf R)$ is localized. With the gauge-locked `filukk` it is:

![real-space decay & Koster-Slater truncation](../assets/tmat_decay_locality.png)
*On-site $\|\tilde V^W(\mathbf R,\mathbf R)\|$ decays $\sim$370× (V$_S$) / 160× (O$_S$) over $\lesssim$5 Å (envelope $\lambda\approx1.5$/3.0 Å); the Koster–Slater truncation captures $\geq$99 % already at $R_{\rm cut}{=}0$ and 100 % by $R_{\rm cut}{=}2$ — so $R_{\rm cut}{=}4$ is amply converged.*

![electron-index gauge check](../assets/tmat_decay_gauge.png)
*Gauge diagnostic: the electron-index decay $\|M^W(R_e;q)\|$ now falls for **every** $q$ (48–79× over $\sim$5 Å). A per-$\mathbf k$ gauge mismatch would leave $q\neq0$ flat (the failure mode of the initial stale `filukk`, where $R_{\rm cut}{=}4$ captured only 55 %). All-$q$ decay confirms $U(\mathbf k)$ and $M$ share one smooth gauge.*

## 7. Summary

| Quantity | V$_S$ (vacancy) | O$_S$ (substitution) |
|---|---|---|
| Downfolded in-gap level (vs Anvil) | $+1.21$ eV ✓ ($+1.205$) | $+0.73$ eV ✓ ($+0.731$) |
| $\|T_{PP}\|$ at K (VBM) | $\approx0.004$ Ry | $\approx0.080$ Ry (**$\sim$19×**) |
| $\|T_{PP}\|$ max on path | $\approx0.035$ Ry | $\approx0.080$ Ry |
| VBM broadening $\|{\rm Im}\,n_d\Sigma\|$ max ($n_d=2.78\%$) | $\sim$9.5 meV | $\sim$13 meV |
| VBM level shift Re $n_d\Sigma$ | $-1$ to $+13$ meV | $-16$ to $-30$ meV |
| Rest-space dressing | non-trivial | mild |
| $\tilde V^W$ locality ($\lambda$; $R_{\rm cut}{=}4$ capture) | 1.5 Å; 100 % | 3.0 Å; 100 % |

**Headline:** at the valence-band edge (K) the isovalent **O$_S$ scatters $\sim$19× harder than the V$_S$ vacancy** and shifts the VBM down markedly more, while the **rest-space dressing matters more for V$_S$** — the beyond-Born physics the explicit-summation $T$-matrix resolves. All quantities use the gauge-locked Wannier rotation (validated by the §6 real-space decay). Pipeline and matrix elements are in place to repeat this for WS₂ V$_S$/O$_S$.
