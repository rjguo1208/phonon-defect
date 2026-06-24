# Static vs dynamic rest-space in the downfolded $T$-matrix

### Derivation: the rest dressing $T^R$, its static limit at $\omega_0$, the exact frequency-dependent form, and the validity condition $W_A/\Delta\ll1$

The electron–defect $T$-matrix is built in two layers (active/rest, cRPA-style): the **rest** bands renormalize the bare defect potential into an effective potential $T^R$, then the **active** multiple scattering is resummed exactly. The production pipeline freezes the rest dressing at a single reference energy $\omega_0$ — a **static approximation**. This note derives where that approximation enters, writes the exact (frequency-dependent) alternative, and quantifies the error: it is controlled by $W_A/\Delta$ (active bandwidth over the active–rest gap), so it degrades as the active space grows.

## 1. Two-layer $T$-matrix and the static approximation

With $g$ the bare defect potential ($V$ / the EDI matrix element $M$), $G^R$ the **rest**-space host propagator (bands 18–70), and $G^A$ the **active**-space host propagator (the 11 Wannier bands 7–17):

$$T^R=[\,1-g\,G^R(\omega_0)\,]^{-1}g,\qquad T^A(\omega)=[\,1-T^R\,G^A(\omega)\,]^{-1}T^R.$$

Layer 1 ($T^R$) folds the rest bands into a dressed potential; Layer 2 ($T^A$) resums the active dynamics. The **static approximation** is the freezing of the rest dressing at $\omega_0$ (here $\omega_0=-5.955$ eV, near the VBM): $T^R\equiv T^R(\omega_0)$, computed **once**, while $G^A(\omega)$ is kept fully dynamic.

## 2. Bridge to the Feshbach self-energy

Projecting $T^R$ onto the active block $P$ (block algebra with $G^R=Q(\omega-E_Q)^{-1}Q$):

$$P\,T^R P=g_{PP}+g_{PQ}\,G^R\,[1-g_{QQ}G^R]^{-1}g_{QP}
=g_{PP}+g_{PQ}(\omega_0-E_Q-g_{QQ})^{-1}g_{QP}=W_{PP}+\Sigma_P(\omega_0)\equiv\tilde V.$$

So the rest $T$-matrix (resummed to **all orders** in the rest self-block $g_{QQ}$) is *identically* the Feshbach downfolded effective potential $\tilde V=W_{PP}+\Sigma_P$, with

$$\Sigma_P(\omega)=W_{PQ}\,(\omega-H_{QQ})^{-1}W_{QP}=B^\dagger\,\mathrm{diag}\!\Big(\tfrac{1}{\omega-\lambda_r}\Big)B,\qquad H_{QQ}=U\,\mathrm{diag}(\lambda_r)\,U^\dagger,\ \ B=U^\dagger W_{QP}.$$

The two languages — rest $T$-matrix and Feshbach self-energy — are the same object. (Here $W=M\,\mathrm{Ry}/N_{\mathbf k}$, the physically normalized matrix element.)

## 3. The two-layer split is exact — iff $T^R$ runs with $\omega$

Because the host $H_0$ is diagonal in the band basis, $G_0=G^A\oplus G^R$ (no active–rest cross terms). The full $T$-matrix $T(\omega)=[1-g\,G_0(\omega)]^{-1}g$ factorizes exactly. Using

$$1-g(G^A+G^R)=(1-gG^R)\,\big[\,1-(1-gG^R)^{-1}g\,G^A\,\big],\qquad (1-gG^R)^{-1}g=T^R,$$

one obtains

$$\boxed{\;T(\omega)=[\,1-T^R(\omega)\,G^A(\omega)\,]^{-1}T^R(\omega),\qquad T^R(\omega)=[\,1-g\,G^R(\omega)\,]^{-1}g\;}$$

This is an **exact identity** (resum rest, then active = resum $G_0$ at once) — *provided $T^R$ is evaluated at the same running $\omega$ as $G^A$*. Therefore the static replacement $T^R(\omega)\!\to\!T^R(\omega_0)$ is the **only** frequency approximation in Layer 1; Layer 2 is exact. Remove it (and take enough rest bands) and the two-layer scheme returns the exact full $T$-matrix.

## 4. The exact (dynamic) rest dressing

Removing the static approximation means re-evaluating the rest dressing at every probe frequency:

$$\boxed{\;T^A(\omega)=[\,1-T^R(\omega)\,G^A(\omega)\,]^{-1}T^R(\omega),\qquad
T^R(\omega)=W_{PP}+\Sigma_P(\omega),\ \ \Sigma_P(\omega)=B^\dagger\,\mathrm{diag}\!\Big(\tfrac{1}{\omega-\lambda_r}\Big)B\;}$$

Because the eigendecomposition $(B,\lambda_r)$ of $H_{QQ}$ is computed once, $\Sigma_P(\omega)$ at any $\omega$ costs a single matrix product — so going dynamic is cheap; the only added work per $\omega$ is recomputing $T^R(\omega)$ and its rotation into the defect block before the Koster–Slater inversion.

## 5. Error of the static approximation and validity condition

Differentiating $T^R=[1-gG^R]^{-1}g$ (with $\partial_\omega A^{-1}=-A^{-1}(\partial_\omega A)A^{-1}$, $\partial_\omega A=-g\,\partial_\omega G^R$):

$$\partial_\omega T^R=T^R\,(\partial_\omega G^R)\,T^R,\qquad
\partial_\omega G^R(\omega)=-\sum_r\frac{|r\rangle\langle r|}{(\omega-\lambda_r)^2}\preceq0.$$

Hence the static error to leading order,

$$\delta T^R=T^R(\omega)-T^R(\omega_0)\approx(\omega-\omega_0)\,T^R(\partial_\omega G^R)T^R,\qquad \|\partial_\omega G^R\|\sim\frac{1}{\Delta^2},$$

so the relative error scales as $(\omega-\omega_0)\,\|T^RG^R\|/\Delta\sim W_A/\Delta$:

$$\boxed{\ \text{static approximation valid}\iff W_A/\Delta\ll1\ }$$

with $W_A$ the active frequency window of interest and $\Delta$ the gap from the active window to the rest bands (18+). Enlarging the active space — or pushing the $P/Q$ boundary up — increases $W_A$ and shrinks $\Delta$, so the static dressing degrades. This is the precise statement of the intuition that a large active space breaks the static approximation.

**Equivalently, in self-energy / $Z$-factor form.** With $\beta_n=-\langle\psi_n|\Sigma_P'(\omega)|\psi_n\rangle=\sum_r|(B\psi_n)_r|^2/(\omega-\lambda_r)^2\ge0$ and $Z_n=(1+\beta_n)^{-1}\le1$, the static level error is

$$\omega_n^{\rm stat}-\omega_n\approx(1-Z_n)\,[\lambda_n(\omega_0)-\omega_0]\approx\beta_n\,(\omega_n-\omega_0),$$

i.e. it vanishes for a level sitting at $\omega_0$ and grows linearly with its distance from $\omega_0$ (worst case $\sim W_A$), weighted by the self-energy slope $\beta_n\sim\|W_{PQ}\|^2/\Delta^2$.

## 6. Where the static error shows up (and where it doesn't)

The static approximation acts only on $T^R$ (Layer 1); the bound-state poles are set mainly by the dynamic $G^A(\omega)$ (Layer 2, exact). So:

- **Pole positions are weakly affected.** For MoS$_2$ V$_S$ the $e$ level sits $\sim1.21$ eV from $\omega_0$, yet static ($+1.209$) vs self-consistent/dynamic ($+1.192$) differ by only $\sim17$ meV — $Z\approx0.99$ for this 11-band partition.
- **Spectral weight and satellites are more affected.** A static $T^R$ gives the active quasiparticle peaks weight $1$ (over-weighted: misses the $Z<1$ renormalization), and omits the rest/satellite structure, so $\int A(\mathbf k,\omega)\,d\omega$ is not exactly satisfied. The spectral function exposes the static approximation more than any single pole.
- **Large active space.** Both pole and weight errors grow as $W_A/\Delta$.

| Treatment | Layer-1 dressing | Cost | Error |
|---|---|---|---|
| **Static (production)** | $T^R(\omega_0)$, computed once | 1 dressing | $O(W_A/\Delta)$ |
| Linearized + $Z$ | $T^R(\omega_0)+(\omega-\omega_0)\partial_\omega T^R$ | + one $\partial_\omega T^R$ | $O((W_A/\Delta)^2)$ |
| **Dynamic (exact)** | $T^R(\omega)$ at each $\omega$ | one GEMM per $\omega$ | exact (within rest-band truncation) |

## 7. Proposed comparison

On the $\Gamma$–M–K path, contrast static $T^R(\omega_0)$ vs dynamic $T^R(\omega)$ for: (1) defect-pole position vs $|\omega_n-\omega_0|$; (2) the spectral function $A(\mathbf k,\omega)$ (peak shift, $Z$ weight, satellites, $\int A$ sum rule); (3) a sweep of the active-space size (11 → 15 → 21 bands) to trace the static error against $W_A/\Delta$. Expectation: with the current 11-band, large-$\Delta$ partition the static dressing is good to $\sim$15 meV on the pole but already shows weight deviations; enlarging the active space degrades both as $W_A/\Delta$. *(Derivation only — not yet run.)*
