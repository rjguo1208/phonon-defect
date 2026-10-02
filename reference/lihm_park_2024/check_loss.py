import numpy as np
rng = np.random.default_rng(7)
Nph = 3; w_ph = np.array([33.0, 35.0, 60.0]); U = 2.0
g0 = 6.0; g_ph = np.array([0.4, 0.3, 0.8]); w0, Om0 = 40.0, 30.0
glr = 8.0 * (rng.normal(size=Nph) + 1j * rng.normal(size=Nph))
gq0 = Om0 * np.sqrt(U); c = gq0 * glr / U
D0 = lambda w, wa, ga: 1.0 / ((w + 1j * ga) ** 2 - wa ** 2)
N = Nph + 1; Gam = np.diag(np.r_[g0, g_ph])
C = np.diag(np.r_[w0, w_ph] ** 2).astype(complex); C[0, 1:] = c; C[1:, 0] = c.conj()
Dt = lambda w: np.linalg.inv((w * np.eye(N) + 1j * Gam) @ (w * np.eye(N) + 1j * Gam) - C)
D = lambda w: np.linalg.inv(np.diag(1 / D0(w, w_ph, g_ph)) - np.outer(c.conj(), c) * D0(w, w0, g0))
# long-range W_00 (S72-S75): vertices g^lr_q0 = gq0 (plasmon), g^lr_qnu = glr (phonons)
gv = np.r_[gq0, glr]
ws = np.linspace(0.01, 150.0, 30000); dw = ws[1] - ws[0]
Wup = np.array([U + gv @ Dt(w) @ gv.conj() for w in ws])                          # S75
def Wlit(w):
    d0 = D0(w, w0, g0); L = glr + gq0 * d0 * c                                    # S55 with lr vertices
    return U + gq0 * d0 * gq0 + L @ D(w) @ np.conj(L)                             # S73 + S74 literal
Wl = np.array([Wlit(w) for w in ws])
for name, W in [("upfolded S75", Wup), ("S74 literal conj", Wl)]:
    loss = -W.imag / U                                                            # ∝ -Im eps^-1_pl-ph
    print(f"{name:18s}: min(-Im W)/U = {loss.min():+.3e}   integral = {loss.sum()*dw:.4f}   "
          f"peak at {ws[loss.argmax()]:.2f} meV")

# independent RPA truth: W = W_b / (1 - dchi0 W_b), W_b = U + sum_nu |g_nu|^2 D0_nu (undoped, static conj only)
def Wtruth(w):
    Ud = Om0 ** 2 * D0(w, w0, g0)              # U*dchi (interacting, plasmon-pole S3)
    dchi0 = (Ud / U) / (1 + Ud)                # dchi0 = dchi / (1 + U dchi)
    Wb = U + np.sum(np.abs(glr) ** 2 * D0(w, w_ph, g_ph))
    return Wb / (1 - dchi0 * Wb)
Wt = np.array([Wtruth(w) for w in ws])
print(f"max |W_up  - W_RPA|/|W_RPA| = {np.max(np.abs(Wup - Wt) / np.abs(Wt)):.2e}")
print(f"max |W_lit - W_RPA|/|W_RPA| = {np.max(np.abs(Wl - Wt) / np.abs(Wt)):.2e}")
