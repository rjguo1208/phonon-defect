"""Numerical check of the plasmon-phonon 'upfolding' in Lihm & Park, arXiv:2409.07393 (SM Sec. S1-S2).

Toy model: Nph phonons + 1 plasmon pole, random complex vertices, finite linewidths.
Checks
  (1) S14->S15 : phonon block of the upfolded D~ equals the direct Dyson solution  [D0^-1 - Pi]^-1
  (2) S23/Eq.9 : D~(w) from the w-independent 2(N+1) matrix K = [[-iG, 1],[C, -iG]] via
                 biorthogonal eigen-expansion (one diagonalisation, all w)
  (3) S32/S33  : spectrum of K comes in pairs (lam, -lam*), Im lam < 0 (retarded)
  (4) S59      : W = U + sum g_a D~_ab g_b^*  (UNSCREENED static vertices in the upfolded space)
                 equals W_el + W_ph built from the w-dependent SCREENED vertex (S55-S58)
  (5) S49      : multipole with a pole-independent c (as printed) vs the pole-resolved c_p
"""
import numpy as np

rng = np.random.default_rng(7)
Nph = 3
w_ph = np.array([33.0, 35.0, 60.0])          # meV
g_ph = np.array([0.4, 0.3, 0.8])              # phonon linewidths
U = 2.0                                       # screened Coulomb U_q
w0, Om0, g0 = 40.0, 30.0, 6.0                 # plasmon frequency, strength, linewidth
glr = rng.normal(size=Nph) + 1j * rng.normal(size=Nph)    # long-range e-ph vertex g^lr_qnu
gq0 = Om0 * np.sqrt(U)                        # S8: long-range e-plasmon vertex (real)
c = gq0 * glr / U                             # S10 coupling amplitude (= Om0 glr / sqrt(U))

def D0(w, wa, ga):
    return 1.0 / ((w + 1j * ga) ** 2 - wa ** 2)

def Pi(w):                                    # S11
    return np.outer(c.conj(), c) * D0(w, w0, g0)

def D_dyson(w):                               # S12
    return np.linalg.inv(np.diag(1 / D0(w, w_ph, g_ph)) - Pi(w))

N = Nph + 1
Gam = np.diag(np.r_[g0, g_ph])
C = np.zeros((N, N), complex)                 # S18
C[0, 0] = w0 ** 2
C[np.arange(1, N), np.arange(1, N)] = w_ph ** 2
C[0, 1:] = c
C[1:, 0] = c.conj()

def Dt_direct(w):                             # S22
    a = (w * np.eye(N) + 1j * Gam)
    return np.linalg.inv(a @ a - C)

# (2) w-independent 2N x 2N matrix (main text Eq. 9) and its biorthogonal eigen-expansion
K = np.block([[-1j * Gam, np.eye(N)], [C, -1j * Gam]])
lam, R = np.linalg.eig(K)
L = np.linalg.inv(R)                          # rows = left eigenvectors, L R = 1

def Dt_modes(w):                              # (1 0) (w - K)^-1 (0 1)^T
    return (R[:N, :] * (1.0 / (w - lam))) @ L[:, N:]

ws = np.linspace(0.0, 120.0, 241)
e1 = max(np.abs(Dt_direct(w)[1:, 1:] - D_dyson(w)).max() / np.abs(D_dyson(w)).max() for w in ws)
e2 = max(np.abs(Dt_modes(w) - Dt_direct(w)).max() / np.abs(Dt_direct(w)).max() for w in ws)
print(f"(1) phonon block of upfolded D~ vs direct Dyson : max rel err = {e1:.2e}")
print(f"(2) eigen-mode expansion of K vs direct inverse  : max rel err = {e2:.2e}")

# (3) pairing lam <-> -lam*
lam_s = np.sort_complex(lam)
pair_err = min(np.abs(np.sort_complex(-lam.conj()) - lam_s).max(), 1e9)
print(f"(3) Im(lam) max = {lam.imag.max():.3f} (<0 => retarded);  |{{lam}} - {{-lam*}}| = {pair_err:.2e}")
pos = lam[lam.real > 0]
print("    hybrid modes  w~ - i g~ :", ", ".join(f"{z.real:.2f}{z.imag:+.2f}i" for z in np.sort_complex(pos)))

# (4) electron-hybrid coupling through UNSCREENED static vertices
g_sr = rng.normal(size=Nph) + 1j * rng.normal(size=Nph)
ov = 0.9 * np.exp(0.3j)                       # <u_mk+q|u_nk>
g_mn0 = gq0 * ov                              # S53
g_mnnu = g_sr + glr * ov                      # S54 (unscreened e-ph vertex)
g_vec = np.r_[g_mn0, g_mnnu]
e4_lit = e4_ana = 0.0
for w in ws:
    W_up = U + g_vec @ Dt_direct(w) @ g_vec.conj()                     # S59
    d0 = D0(w, w0, g0)
    gscr = g_mnnu + g_mn0 * d0 * c                                     # S55
    W_el = U + g_mn0 * d0 * np.conj(g_mn0)                             # S57
    W_lit = W_el + gscr @ D_dyson(w) @ gscr.conj()                     # S58 read literally: conj also hits eps^-1(w)
    W_ana = W_el + gscr @ D_dyson(w) @ (g_mnnu.conj() + np.conj(g_mn0) * d0 * c.conj())  # conj static parts only
    e4_lit = max(e4_lit, abs(W_up - W_lit) / abs(W_up))
    e4_ana = max(e4_ana, abs(W_up - W_ana) / abs(W_up))
print(f"(4) S59 vs S58, '*' conjugating static vertex only : max rel err = {e4_ana:.2e}")
print(f"    S59 vs S58, '*' read literally (conj eps^-1(w)) : max rel err = {e4_lit:.2e}  <- wrong for finite plasmon linewidth")

# (5) multipole: two poles with different strengths
wp = np.array([25.0, 70.0]); Omp = np.array([18.0, 40.0]); gp = np.array([4.0, 9.0])
def Pi_true(w):                                # from S47: dchi = (1/U) sum_p Om_p^2 D0_p
    s = np.sum(Omp ** 2 * D0(w, wp, gp)) / U
    return np.outer(glr.conj(), glr) * s
def build(cp_rows):
    Np = len(wp); M = Np + Nph
    Cm = np.zeros((M, M), complex)
    Cm[np.arange(Np), np.arange(Np)] = wp ** 2
    Cm[np.arange(Np, M), np.arange(Np, M)] = w_ph ** 2
    Cm[:Np, Np:] = cp_rows; Cm[Np:, :Np] = cp_rows.conj().T
    G = np.diag(np.r_[gp, g_ph])
    return lambda w: np.linalg.inv((w * np.eye(M) + 1j * G) @ (w * np.eye(M) + 1j * G) - Cm)[Np:, Np:]
c_pole = Omp[:, None] * glr[None, :] / np.sqrt(U)          # c_{p nu} = Om_p g^lr_nu / sqrt(U)
c_flat = np.tile(Om0 * glr / np.sqrt(U), (len(wp), 1))     # as printed in S49-S51 (no p index)
D_true = lambda w: np.linalg.inv(np.diag(1 / D0(w, w_ph, g_ph)) - Pi_true(w))
for name, cp in [("pole-resolved c_p", c_pole), ("pole-independent c (as printed)", c_flat)]:
    f = build(cp)
    err = max(np.abs(f(w) - D_true(w)).max() / np.abs(D_true(w)).max() for w in ws)
    print(f"(5) multipole, {name:32s}: max rel err vs direct Dyson = {err:.2e}")
