# Fan–Migdal self-energy and the phonon spectral function in a defective crystal

This note contains two linked derivations. **Part I** (§1–§3) derives the Fan–Migdal (FM) self-energy — the lowest non-vanishing term of the electron–phonon (e–ph) interaction — and shows explicitly that it is built from the *bare* electron Green's function $G^0$ and the *bare* phonon propagator $D^0$. **Part II** (§4–§8) derives how the *phonon* spectral function is computed when the crystal contains a dilute concentration of point defects: an exact single-defect $T$-matrix in a small defect-localized subspace, a configurational average that restores crystal momentum, the resulting $B_{\mathbf q\nu}(\omega)$, and the electron-mediated (bubble) contribution with its double-counting subtraction. Units: $\hbar=k_B=1$. $N_p$ is the number of primitive cells in the Born–von-Kármán supercell; fermionic Matsubara frequencies $i\omega_j=(2j{+}1)\pi/\beta$, bosonic $iq_l=2l\pi/\beta$, $\beta=1/T$.

## 1. Setup: Hamiltonian and bare propagators

The standard ab-initio e–ph Hamiltonian, with electrons in Kohn–Sham (KS) bands $\varepsilon_{n\mathbf k}$, harmonic phonons $\omega_{\mathbf q\nu}$, and the *linear* coupling:

$$
\hat H=\underbrace{\sum_{n\mathbf k}\varepsilon_{n\mathbf k}\,\hat c^\dagger_{n\mathbf k}\hat c_{n\mathbf k}}_{\hat H_e}
+\underbrace{\sum_{\mathbf q\nu}\omega_{\mathbf q\nu}\Big(\hat a^\dagger_{\mathbf q\nu}\hat a_{\mathbf q\nu}+\tfrac12\Big)}_{\hat H_{\rm ph}}
+\underbrace{\frac{1}{\sqrt{N_p}}\sum_{\mathbf k\mathbf q}\sum_{mn\nu}g_{mn\nu}(\mathbf k,\mathbf q)\,
\hat c^\dagger_{m\mathbf k+\mathbf q}\hat c_{n\mathbf k}\,\hat A_{\mathbf q\nu}}_{\hat H_{ep}},
\qquad
\hat A_{\mathbf q\nu}\equiv\hat a_{\mathbf q\nu}+\hat a^\dagger_{-\mathbf q\nu},
$$

where the e–ph vertex is the first-order variation of the self-consistent KS potential along the phonon eigendisplacement, with the zero-point amplitude pulled in:

$$
g_{mn\nu}(\mathbf k,\mathbf q)=\Big(\frac{1}{2\omega_{\mathbf q\nu}}\Big)^{1/2}
\big\langle u_{m\mathbf k+\mathbf q}\big|\Delta_{\mathbf q\nu}v^{\rm KS}\big|u_{n\mathbf k}\big\rangle ,
\qquad
g_{nm\nu}(\mathbf k{+}\mathbf q,-\mathbf q)=g^{*}_{mn\nu}(\mathbf k,\mathbf q).
$$

The *bare* (non-interacting) finite-temperature propagators are defined with respect to $\hat H_0=\hat H_e+\hat H_{\rm ph}$:

$$
G^0_{n\mathbf k}(\tau)=-\big\langle T_\tau\,\hat c_{n\mathbf k}(\tau)\,\hat c^\dagger_{n\mathbf k}(0)\big\rangle_0
\;\;\Longrightarrow\;\;
G^0_{n\mathbf k}(i\omega_j)=\frac{1}{i\omega_j-\xi_{n\mathbf k}},\qquad \xi_{n\mathbf k}\equiv\varepsilon_{n\mathbf k}-\mu ,
$$

$$
D^0_{\mathbf q\nu}(\tau)=-\big\langle T_\tau\,\hat A_{\mathbf q\nu}(\tau)\,\hat A^\dagger_{\mathbf q\nu}(0)\big\rangle_0 ,
\qquad \hat A^\dagger_{\mathbf q\nu}=\hat A_{-\mathbf q\nu}.
$$

The phonon propagator is worth deriving once. In the Matsubara Heisenberg picture $\hat a_{\mathbf q\nu}(\tau)=\hat a_{\mathbf q\nu}e^{-\omega_{\mathbf q\nu}\tau}$, $\hat a^\dagger_{\mathbf q\nu}(\tau)=\hat a^\dagger_{\mathbf q\nu}e^{+\omega_{\mathbf q\nu}\tau}$, and with $n_{\mathbf q\nu}=(e^{\beta\omega_{\mathbf q\nu}}-1)^{-1}$, for $0<\tau<\beta$:

$$
D^0_{\mathbf q\nu}(\tau)=-\big[(n_{\mathbf q\nu}+1)\,e^{-\omega_{\mathbf q\nu}\tau}+n_{\mathbf q\nu}\,e^{+\omega_{\mathbf q\nu}\tau}\big].
$$

Fourier transforming, $D^0(iq_l)=\int_0^\beta d\tau\,e^{iq_l\tau}D^0(\tau)$, and using $e^{iq_l\beta}=1$ together with $(n+1)(e^{-\beta\omega}-1)=-1$ and $n(e^{\beta\omega}-1)=1$:

$$
D^0_{\mathbf q\nu}(iq_l)=\frac{1}{iq_l-\omega_{\mathbf q\nu}}-\frac{1}{iq_l+\omega_{\mathbf q\nu}}
=\frac{2\omega_{\mathbf q\nu}}{(iq_l)^2-\omega_{\mathbf q\nu}^2}.
$$

These two objects — and nothing else — will appear inside the lowest-order self-energy.

## 2. Perturbation expansion: the lowest non-vanishing term

Expand the full electron Green's function in powers of $\hat H_{ep}$ with the imaginary-time $S$-matrix:

$$
G_{n\mathbf k}(\tau)=-\frac{\Big\langle T_\tau\,\hat c_{n\mathbf k}(\tau)\,\hat c^\dagger_{n\mathbf k}(0)\,
\exp\big[-\int_0^\beta d\tau'\,\hat H_{ep}(\tau')\big]\Big\rangle_0}
{\Big\langle T_\tau\exp\big[-\int_0^\beta d\tau'\,\hat H_{ep}(\tau')\big]\Big\rangle_0}.
$$

**First order vanishes.** The $O(g)$ term contains a single phonon operator, and $\langle\hat A_{\mathbf q\nu}\rangle_0=0$: in thermal equilibrium with the ions at their relaxed positions the lattice carries no static distortion (any nonzero $\langle\hat A\rangle$ would simply redefine the reference geometry and is absorbed there). The first non-vanishing term is therefore *second* order — one phonon emitted and reabsorbed.

**Second order: Wick contractions.** At $O(g^2)$, Wick's theorem factorizes the expectation value into fermion and boson contractions. The boson contraction ties the two vertices into a single phonon line; momentum conservation forces the two vertices to carry $\mathbf q$ and $-\mathbf q$:

$$
\big\langle T_\tau\,\hat A_{\mathbf q\nu}(\tau_1)\,\hat A_{\mathbf q'\nu'}(\tau_2)\big\rangle_0
=-\,\delta_{\mathbf q',-\mathbf q}\,\delta_{\nu\nu'}\;D^0_{\mathbf q\nu}(\tau_1-\tau_2).
$$

The fermion contractions give two topologies: (a) the **connected "rainbow" diagram**, where a bare electron line $G^0_{m\mathbf k+\mathbf q}$ propagates between the two vertices under the phonon arc; (b) **tadpole / disconnected pieces** — a closed fermion loop attached to one vertex by a $\mathbf q=0$ phonon line. The disconnected pieces cancel against the denominator (linked-cluster theorem), and the $\mathbf q=0$ tadpole, $\propto\sum_{n\mathbf k}g_{nn\nu}(\mathbf k,0)f_{n\mathbf k}$, is the first-order force on the lattice — it vanishes by the same equilibrium-geometry condition that killed the first-order term. The $1/2!$ of the exponential is cancelled by the two equivalent vertex orderings. Collecting signs from the contractions, what survives is

$$
G^{(2)}_{n\mathbf k}(\tau)=\int_0^\beta\!\!d\tau_1\!\int_0^\beta\!\!d\tau_2\;
G^0_{n\mathbf k}(\tau-\tau_1)\,\Sigma^{\rm FM}_{n\mathbf k}(\tau_1-\tau_2)\,G^0_{n\mathbf k}(\tau_2),
$$

$$
\Sigma^{\rm FM}_{n\mathbf k}(\tau)=-\frac{1}{N_p}\sum_{m\nu}\sum_{\mathbf q}
\big|g_{mn\nu}(\mathbf k,\mathbf q)\big|^2\;G^0_{m\mathbf k+\mathbf q}(\tau)\;D^0_{\mathbf q\nu}(\tau),
$$

i.e. in frequency space (the product in $\tau$ becomes a convolution):

$$
\boxed{\;\;
\Sigma^{\rm FM}_{n\mathbf k}(i\omega_j)\;=\;-\frac{1}{\beta N_p}\sum_{m\nu}\sum_{\mathbf q}\sum_{iq_l}
\big|g_{mn\nu}(\mathbf k,\mathbf q)\big|^2\;
G^0_{m\mathbf k+\mathbf q}(i\omega_j+iq_l)\;D^0_{\mathbf q\nu}(iq_l).
\;\;}
$$

Both internal lines are **bare**: $G^0$ for the electron, $D^0$ for the phonon — exactly the statement to be proven. Resumming the one-particle-irreducible part to all orders gives the Dyson equation $G^{-1}_{n\mathbf k}(i\omega_j)=i\omega_j-\xi_{n\mathbf k}-\Sigma_{n\mathbf k}(i\omega_j)$, of which the boxed expression is the leading skeleton.

**Why this is the right "lowest order" (Migdal).** Corrections that dress the e–ph *vertex* are smaller by $O(\omega_{\rm ph}/E_F)\sim O(\sqrt{m/M})$ (Migdal's theorem), so the bare-vertex rainbow is the dominant skeleton. Self-consistently dressing the internal $G$ or $D$ modifies the result only at the same order as the neglected vertex corrections; the one-shot $G^0D^0$ evaluation above is therefore the standard, internally consistent lowest-order theory.

**Remark (Debye–Waller).** At the same order in the ionic displacement, $O(u^2)$, there is a second diagram: the *Debye–Waller* self-energy, which is first order in the *quadratic* e–ph vertex $\partial^2 v^{\rm KS}/\partial u^2$. It is static and real (it shifts band energies but carries no linewidth). The full $O(u^2)$ theory (Allen–Heine–Cardona) is $\Sigma=\Sigma^{\rm FM}+\Sigma^{\rm DW}$; the FM term is the only one built from the linear coupling and the only dynamical one. We focus on $\Sigma^{\rm FM}$.

## 3. The Matsubara sum, analytic continuation, and the physical content

Evaluate the bosonic frequency sum exactly. Split $D^0$ into partial fractions, $D^0(iq_l)=(iq_l-\omega_{\mathbf q\nu})^{-1}-(iq_l+\omega_{\mathbf q\nu})^{-1}$, and use the standard contour result for bosonic sums (poles of $n_B(z)=(e^{\beta z}-1)^{-1}$ at $z=iq_l$ with residue $1/\beta$):

$$
\frac1\beta\sum_{l}F(iq_l)=-\sum_{z_p}n_B(z_p)\,\mathrm{Res}_{z=z_p}F(z),
\qquad
F(z)=\frac{1}{z-a}\cdot\frac{1}{z+i\omega_j-\xi}.
$$

$F$ has simple poles at $z=a$ (residue $1/(a+i\omega_j-\xi)$) and at $z=\xi-i\omega_j$ (residue $1/(\xi-i\omega_j-a)$). Because $i\omega_j$ is fermionic, $e^{-i\beta\omega_j}=-1$, so $n_B(\xi-i\omega_j)=(e^{\beta\xi}e^{-i\beta\omega_j}-1)^{-1}=-n_F(\xi)$ — a boson distribution evaluated one fermionic frequency away turns into (minus) a Fermi function. Hence

$$
-\frac1\beta\sum_{l}\frac{1}{iq_l-a}\cdot\frac{1}{iq_l+i\omega_j-\xi}
=\frac{n_B(a)+n_F(\xi)}{i\omega_j+a-\xi}.
$$

Apply this with $a=+\omega_{\mathbf q\nu}$ and $a=-\omega_{\mathbf q\nu}$, using $n_B(-\omega)=-[1+n_B(\omega)]$:

$$
\Sigma^{\rm FM}_{n\mathbf k}(i\omega_j)=\frac{1}{N_p}\sum_{m\nu}\sum_{\mathbf q}
\big|g_{mn\nu}(\mathbf k,\mathbf q)\big|^2
\left[\frac{n_{\mathbf q\nu}+f_{m\mathbf k+\mathbf q}}{i\omega_j-\xi_{m\mathbf k+\mathbf q}+\omega_{\mathbf q\nu}}
+\frac{n_{\mathbf q\nu}+1-f_{m\mathbf k+\mathbf q}}{i\omega_j-\xi_{m\mathbf k+\mathbf q}-\omega_{\mathbf q\nu}}\right],
$$

with $f$ and $n$ the Fermi–Dirac and Bose–Einstein occupations at temperature $T$. Analytic continuation $i\omega_j\to\omega+i\eta$ gives the retarded self-energy; its imaginary part is

$$
\mathrm{Im}\,\Sigma^{\rm FM}_{n\mathbf k}(\omega)=-\frac{\pi}{N_p}\sum_{m\nu\mathbf q}|g_{mn\nu}(\mathbf k,\mathbf q)|^2
\Big[(n_{\mathbf q\nu}+f_{m\mathbf k+\mathbf q})\,\delta(\omega-\xi_{m\mathbf k+\mathbf q}+\omega_{\mathbf q\nu})
+(n_{\mathbf q\nu}+1-f_{m\mathbf k+\mathbf q})\,\delta(\omega-\xi_{m\mathbf k+\mathbf q}-\omega_{\mathbf q\nu})\Big].
$$

Read on the quasiparticle shell $\omega=\xi_{n\mathbf k}$ this is exactly Fermi's golden rule: the first term is phonon **absorption** ($\propto n$), the second **emission** ($\propto n+1$, spontaneous + stimulated), with Pauli blocking carried by the $f$ factors. The quasiparticle linewidth (inverse lifetime) is $\Gamma_{n\mathbf k}=-2\,\mathrm{Im}\,\Sigma^{\rm FM}_{n\mathbf k}(\xi_{n\mathbf k})=1/\tau_{n\mathbf k}$, and the measurable electron spectral function follows from the Dyson-resummed $G$:

$$
A_{n\mathbf k}(\omega)=-\frac1\pi\,\mathrm{Im}\,
\frac{1}{\omega-\xi_{n\mathbf k}-\Sigma^{\rm FM}_{n\mathbf k}(\omega)} .
$$

**What "bare" means in practice.** In ab-initio calculations $G^0$ is built from KS eigenvalues and $D^0$, $g$ from DFPT. The DFPT phonon is *not* the unscreened ion: it already contains the static electronic screening of the ionic motion. Consistency therefore demands that one must **not** additionally dress $D^0$ with the static electronic polarization bubble — that screening is already inside $\omega_{\mathbf q\nu}$ and would be double counted. Only the *non-adiabatic remainder* of the bubble may be added (see §7). This is the precise sense in which the FM self-energy "uses the bare electron and bare phonon Green's functions".

## 4. Part II — Phonons with defects: the exact single-defect T-matrix

Now turn the question around: how does a **defect** modify the **phonon** propagator and spectral function? Work in the harmonic lattice. The host crystal has masses $M_0$ and force constants $\Phi^0_{l\kappa\alpha,l'\kappa'\beta}$; its mass-scaled dynamical operator and normal modes are

$$
\mathcal D_0=M_0^{-1/2}\,\Phi^0\,M_0^{-1/2},
\qquad
\mathcal D_0\,|\mathbf q\nu\rangle=\omega_{\mathbf q\nu}^2\,|\mathbf q\nu\rangle,
\qquad
\langle l\kappa\alpha|\mathbf q\nu\rangle=\frac{e_{\nu\alpha}(\kappa;\mathbf q)}{\sqrt{N_p}}\,e^{i\mathbf q\cdot\mathbf R_l}.
$$

It is most convenient to work with the **resolvent in $z=(\omega+i0^+)^2$**:

$$
g(z)=\big(z-\mathcal D_0\big)^{-1},
\qquad
g_{\mathbf q\nu}(z)=\frac{1}{z-\omega_{\mathbf q\nu}^2}.
$$

A point defect changes both masses and force constants: $M=M_0+\Delta M$, $\Phi=\Phi^0+\Delta\Phi$ (with $\Delta\Phi$ obeying the acoustic sum rule of the *defective* lattice). Scaling the equations of motion $\,\omega^2 M u=\Phi u\,$ by the **host** masses, $\tilde u=M_0^{1/2}u$:

$$
\big[\,z\,(\mathbb 1+\varepsilon)-\mathcal D_0-\Delta\mathcal D\,\big]\,\tilde u=0,
\qquad
\varepsilon=M_0^{-1/2}\Delta M\,M_0^{-1/2},
\quad
\Delta\mathcal D=M_0^{-1/2}\Delta\Phi\,M_0^{-1/2},
$$

so the defective resolvent is $G(z)=[\,z-\mathcal D_0-V(z)\,]^{-1}$ with a **frequency-dependent** perturbation (Lifshitz):

$$
V(z)=\Delta\mathcal D-z\,\varepsilon .
$$

The mass defect enters multiplied by $z$; the force-constant change is static. Crucially, $V(z)$ has support only on the defect site and its few shells of neighbours ($\Delta\Phi$ is short ranged for a neutral defect in a non-polar host; charged defects in polar hosts add a long-range dipole tail that must be treated separately). Dyson's equation and its exact resummation,

$$
G=g+g\,V\,G
\;\;\Longleftrightarrow\;\;
G=g+g\,T\,g,
\qquad
T(z)=\big[\,\mathbb 1-V(z)\,g(z)\,\big]^{-1}V(z),
$$

reduce the problem to inverting $[\mathbb 1-Vg]$ **inside the small defect subspace only** (tens to a few hundred ionic degrees of freedom). This is precisely the phonon analogue of the Koster–Slater electron–defect Green's-function method: $\mathcal D_0$ plays the role of $H_0$, $V(z)$ the role of $\Delta V$, with the one new feature that $V$ depends linearly on $z$ through the mass defect. (For a vacancy, $\Delta M=-M_{\kappa_0}$ and $\Delta\Phi$ decouples the vacant site exactly; its three free degrees of freedom appear as spurious $\omega=0$ solutions of the decoupled block and must be projected out of $T$.)

## 5. Many defects: configurational average and the phonon self-energy

A defect at cell $\mathbf R_l$ carries phases $\langle\mathbf q\nu|T_l|\mathbf q'\nu'\rangle=e^{i(\mathbf q'-\mathbf q)\cdot\mathbf R_l}\,\langle\mathbf q\nu|T_0|\mathbf q'\nu'\rangle$: a single defect breaks momentum conservation. For $N_d$ defects at random positions, averaging over configurations restores it, $\overline{e^{i(\mathbf q'-\mathbf q)\cdot\mathbf R_l}}=\delta_{\mathbf q\mathbf q'}$. In the **dilute, non-overlapping limit** (defect concentration $c=N_d/N_p\ll1$, mean spacing large compared to the range of $V$ and of the single-defect scattered wave) the averaged self-energy is $N_d$ independent single-defect T-matrices — the average-T-matrix approximation (ATA):

$$
\overline G_{\mathbf q\nu}(z)=\frac{1}{z-\omega_{\mathbf q\nu}^2-\pi_{\mathbf q\nu}(z)},
\qquad
\pi_{\mathbf q\nu}(z)=N_d\,\langle\mathbf q\nu|T_0(z)|\mathbf q\nu\rangle\equiv c\;t_{\mathbf q\nu}(z),
$$

where $t_{\mathbf q\nu}\equiv N_p\langle\mathbf q\nu|T_0|\mathbf q\nu\rangle$ is intensive (the $1/N_p$ of the extended-mode normalization cancels the explicit $N_p$). Corrections are $O(c^2)$ (defect pairs) and are systematically resummed, if needed, by the CPA: replace $g\to\overline G$ self-consistently inside $T$. Everything below uses the ATA.

## 6. The phonon spectral function

With the convention $g_{\mathbf q\nu}=(z-\omega_{\mathbf q\nu}^2)^{-1}$ the properly normalized phonon spectral function is

$$
\boxed{\;\;
B_{\mathbf q\nu}(\omega)=-\frac{2\omega}{\pi}\,
\mathrm{Im}\;\frac{1}{(\omega+i0^+)^2-\omega_{\mathbf q\nu}^2-\pi_{\mathbf q\nu}\big((\omega+i0^+)^2\big)}\;,
\qquad
\int_0^\infty\! d\omega\,B_{\mathbf q\nu}(\omega)\simeq1 .
\;\;}
$$

*Normalization check*: for $\pi=0$, $\mathrm{Im}\,[(\omega+i0)^2-\omega_0^2]^{-1}=-\pi\,\delta(\omega^2-\omega_0^2)$ and $-\tfrac{2\omega}\pi\cdot(-\pi)\tfrac{\delta(\omega-\omega_0)}{2\omega_0}=\delta(\omega-\omega_0)$ for $\omega>0$. ✓

**Quasiparticle regime.** When $|\pi|\ll\omega_{\mathbf q\nu}^2$ the line is Lorentzian:

$$
B_{\mathbf q\nu}(\omega)\approx\frac{1}{\pi}\,\frac{\Gamma_{\mathbf q\nu}}{(\omega-\Omega_{\mathbf q\nu})^2+\Gamma_{\mathbf q\nu}^2},
\qquad
\Omega_{\mathbf q\nu}^2=\omega_{\mathbf q\nu}^2+\mathrm{Re}\,\pi_{\mathbf q\nu},
\qquad
\Gamma_{\mathbf q\nu}=-\frac{\mathrm{Im}\,\pi_{\mathbf q\nu}(\Omega^2_{\mathbf q\nu})}{2\,\Omega_{\mathbf q\nu}} ,
$$

so the defect-limited phonon lifetime entering the thermal-conductivity BTE is $1/\tau_{\mathbf q\nu}=2\Gamma_{\mathbf q\nu}$.

**Limit 1 — Born / isotope scattering (consistency check).** Take pure mass disorder, $\Delta\Phi=0$, $V=-z\varepsilon$ on one site $\kappa_0$, and the Born limit $T\to V$. Then $\mathrm{Im}\,\pi_{\mathbf q\nu}=c\,\omega^4\varepsilon^2\,\frac{1}{N_p}\sum_{\mathbf q'\nu'}|\mathbf e^*_{\mathbf q'\nu'}(\kappa_0)\!\cdot\!\mathbf e_{\mathbf q\nu}(\kappa_0)|^2\,\mathrm{Im}\,g_{\mathbf q'\nu'}$, and with $\mathrm{Im}\,g=-\pi\delta(\omega^2-\omega'^2)=-\pi\delta(\omega-\omega')/2\omega$ one finds

$$
\frac{1}{\tau^{\rm iso}_{\mathbf q\nu}}=2\Gamma_{\mathbf q\nu}
=\frac{\pi}{2N_p}\,\omega_{\mathbf q\nu}^2\sum_{\mathbf q'\nu'}\delta(\omega_{\mathbf q\nu}-\omega_{\mathbf q'\nu'})
\sum_{\kappa}g_2(\kappa)\,\big|\mathbf e^{*}_{\mathbf q'\nu'}(\kappa)\cdot\mathbf e_{\mathbf q\nu}(\kappa)\big|^2 ,
\qquad
g_2(\kappa)=\sum_i f_i\Big(1-\frac{m_i(\kappa)}{\bar m(\kappa)}\Big)^{\!2},
$$

which is exactly **Tamura's isotope-scattering formula** — fixing all prefactors of the general result.

**Limit 2 — resonant and localized modes (where Born fails).** The full T-matrix has structure the Born term cannot produce: zeros of $\det\big[\mathbb 1-\mathrm{Re}\,V(z)g(z)\big]$ *inside* the host phonon bands give **quasi-localized resonant modes** — $\pi_{\mathbf q\nu}$ acquires a sharp resonance at $\omega_{\rm res}$, and $B_{\mathbf q\nu}(\omega)$ develops an avoided-crossing/level-repulsion structure pinned at $\omega_{\rm res}$ across all $\mathbf q$ (the flat "defect band" seen in neutron/IXS data and in defective-supercell unfolding). Heavy substitutionals or weakened bonds push $\omega_{\rm res}$ deep into the acoustic branches — these resonances dominate the low-$T$ thermal resistance. Solutions *above* the host spectrum (light defects / stiffened bonds) are true **localized modes**: poles of $T$ at $\omega_{\rm loc}$, appearing in $B_{\mathbf q\nu}$ with weight $O(1/N_p)$ but as sharp peaks in the local spectral function at the defect. The total spectral redistribution obeys the Krein/Friedel trace formula $\Delta\rho(z)=-\frac1\pi\,\mathrm{Im}\,\frac{d}{dz}\ln\det\big[\mathbb 1-V(z)g(z)\big]$, which conserves the total mode count.

## 7. Electronic contribution to the phonon self-energy (and no double counting)

The phonon-side counterpart of the FM diagram is the e–ph **bubble** — the same two vertices, read as a phonon self-energy:

$$
\Pi^{ep}_{\mathbf q\nu}(\omega)=\frac{2}{N_p}\sum_{mn\mathbf k}\big|g_{mn\nu}(\mathbf k,\mathbf q)\big|^2\,
\frac{f_{n\mathbf k}-f_{m\mathbf k+\mathbf q}}{\omega+\varepsilon_{n\mathbf k}-\varepsilon_{m\mathbf k+\mathbf q}+i0^+}
$$

(factor 2 for spin; $\Pi^{ep}$ has units of energy and enters the $z$-resolvent as $\pi^{ep}_{\mathbf q\nu}=2\omega_{\mathbf q\nu}\Pi^{ep}_{\mathbf q\nu}$). As stressed in §3, the DFPT $\omega_{\mathbf q\nu}$ already contain this bubble **statically screened at $\omega=0$**; adding $\Pi^{ep}(\omega)$ wholesale would double count it. The correct prescription is the **non-adiabatic difference**:

$$
\pi^{ep}_{\mathbf q\nu}(\omega)=2\omega_{\mathbf q\nu}\big[\Pi^{ep}_{\mathbf q\nu}(\omega)-\Pi^{ep}_{\mathbf q\nu}(0)\big],
$$

whose imaginary part at $\omega=\Omega_{\mathbf q\nu}$ gives the familiar e–ph phonon linewidth $\gamma_{\mathbf q\nu}$ (Allen). In a **defective** sample this term changes in two ways: (i) defects dope the system — the occupations $f$ and hence the phase space of the bubble become carrier-/defect-concentration dependent (free-carrier Kohn anomalies, phonon softening); (ii) the electron lines inside the bubble should themselves carry the **electron–defect self-energy** — i.e. the bubble is evaluated with defect-broadened spectral functions, $\Gamma^{e\text{-}d}_{n\mathbf k}$ taken from the electron–defect T-matrix/SERTA rates (the quantity the EDI workflow computes), which smears sharp Kohn anomalies on the scale of the electron mean free path. The total phonon self-energy, to leading order in $c$ and $g^2$ and neglecting crossed defect–phonon–electron vertex diagrams (higher order in both),

$$
\pi_{\mathbf q\nu}(\omega)\;=\;
\underbrace{c\,t_{\mathbf q\nu}\big((\omega+i0^+)^2\big)}_{\text{structural defect scattering (§4–5)}}
\;+\;
\underbrace{2\omega_{\mathbf q\nu}\big[\Pi^{ep}_{\mathbf q\nu}(\omega)-\Pi^{ep}_{\mathbf q\nu}(0)\big]}_{\text{non-adiabatic e–ph, defect-dressed electrons}} ,
$$

inserted into the boxed $B_{\mathbf q\nu}(\omega)$ of §6.

## 8. Practical recipe (supercell DFT workflow)

The structure deliberately mirrors the electron–defect (EDI / Koster–Slater) pipeline — supercell difference for the perturbation, host Green's function from the primitive cell, small-block T-matrix:

1. **Host phonons**: DFPT (or finite displacements) on the primitive cell → $\omega_{\mathbf q\nu}$, $\mathbf e_{\mathbf q\nu}$ and Fourier-interpolable force constants $\Phi^0$.
2. **Defect supercell**: relax the defective supercell, compute its force constants $\Phi^{\rm def}$ (finite differences are simplest). Build $\Delta\Phi=\Phi^{\rm def}-\Phi^0$ mapped onto the supercell, and *verify its spatial decay* away from the defect — then truncate to a defect cluster (site + a few neighbour shells). $\Delta M$ lives on the defect site(s); for a vacancy project out the three decoupled degrees of freedom. Re-impose the acoustic sum rule on the truncated $\Delta\Phi$.
3. **Host resolvent on the cluster**: $g_{ab}(z)=\frac{1}{N_p}\sum_{\mathbf q\nu}\frac{\langle a|\mathbf q\nu\rangle\langle\mathbf q\nu|b\rangle\,N_p}{z-\omega_{\mathbf q\nu}^2}$ for cluster orbitals $a,b$, by a dense interpolated $\mathbf q$-sum (tetrahedron method or small $i\eta$) — the analogue of the electronic $g(E)$ in Koster–Slater.
4. **T-matrix**: build $V(z)=\Delta\mathcal D-z\varepsilon$ on the cluster, invert $[\mathbb 1-V(z)g(z)]$ (small dense matrix, one inversion per frequency), project onto host modes → $t_{\mathbf q\nu}(z)$, multiply by the physical defect concentration $c$.
5. **Optional electronic part**: add the non-adiabatic bubble of §7 with defect-dressed electron lines (rates from the electron–defect calculation).
6. **Outputs**: $B_{\mathbf q\nu}(\omega)$ along a $\mathbf q$-path (directly comparable to neutron/IXS and to band-unfolded defective-supercell spectra); mode linewidths $\Gamma_{\mathbf q\nu}\to\tau_{\mathbf q\nu}$ for the phonon BTE (defect-limited lattice thermal conductivity); $\Delta$DOS from the Krein determinant; localized/resonant-mode energies from $\det[\mathbb 1-\mathrm{Re}\,Vg]=0$.

The two boxed results — $\Sigma^{\rm FM}$ from bare $G^0D^0$ (§2–3) and $B_{\mathbf q\nu}(\omega)$ from the defect T-matrix self-energy (§6) — are the deliverables of this note.
