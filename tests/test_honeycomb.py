"""Tests of the honeycomb toy model (run: python -m pytest tests)."""
import os
import sys

import numpy as np
import scipy.sparse.linalg as sla

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from phdef import honeycomb as hc  # noqa: E402


def test_every_term_is_invariant_under_rigid_motions():
    """Acoustic and rotational sum rules: each energy term annihilates all six rigid motions."""
    for atoms, C in hc.terms():
        pts = np.array([hc.position(s, cell) for s, cell in atoms])
        for u in hc.rigid_motions(pts):
            f = [sum(C[a][b] @ u[b] for b in range(len(atoms))) for a in range(len(atoms))]
            assert np.abs(f).max() < 1e-14


def test_bloch_matrix_symmetries():
    rng = np.random.default_rng(0)
    for k in [hc.POINTS["M"], hc.POINTS["K"], *rng.normal(size=(5, 2))]:
        D = hc.bloch(k)
        assert np.allclose(D, D.conj().T, atol=1e-14)
        assert np.abs(D[np.ix_(hc.IN_PLANE, hc.OUT_OF_PLANE)]).max() < 1e-14     # mirror plane
        assert np.allclose(np.linalg.eigvalsh(D), np.linalg.eigvalsh(hc.bloch(k, gauge="cell")), atol=1e-12)
    w2 = np.linalg.eigvalsh(hc.bloch(hc.POINTS["Γ"]))
    assert np.abs(w2[:3]).max() < 1e-13 and w2[3] > 0.5                          # three acoustic zeros


def test_flexural_branch_is_quadratic_and_in_plane_branches_linear():
    ratios = []
    for q in (2e-3, 4e-3, 8e-3):
        D = hc.bloch(np.array([q, 0.0]))
        w_za = np.sqrt(np.linalg.eigvalsh(D[np.ix_(hc.OUT_OF_PLANE, hc.OUT_OF_PLANE)])[0])
        w_ta = np.sqrt(np.linalg.eigvalsh(D[np.ix_(hc.IN_PLANE, hc.IN_PLANE)])[0])
        ratios.append((w_za / q ** 2, w_ta / q))
    za, ta = np.array(ratios).T
    assert np.ptp(za) < 1e-3 * za.mean() and np.ptp(ta) < 1e-4 * ta.mean()


def test_stable_on_a_grid():
    n = 30
    w2 = [np.linalg.eigvalsh(hc.bloch((i * hc.B1 + j * hc.B2) / n)) for i in range(n) for j in range(n)]
    assert np.min(w2) > -1e-12


def test_torus_spectrum_equals_bloch_spectrum_on_the_grid():
    L = 5
    D, _ = hc.supercell(L)
    torus = np.linalg.eigvalsh(D.toarray())
    grid = np.sort(np.concatenate([np.linalg.eigvalsh(hc.bloch((i * hc.B1 + j * hc.B2) / L))
                                   for i in range(L) for j in range(L)]))
    assert np.allclose(torus, grid, atol=1e-12)


def test_substitution_keeps_the_sum_rules():
    L = 5
    D, m = hc.supercell(L, defect=(0, (2, 2), 0.4, 1.3))
    assert abs(D - D.T).max() < 1e-14
    w2 = np.linalg.eigvalsh(D.toarray())
    assert np.abs(w2[:3]).max() < 1e-12 and w2[3] > 1e-3        # still exactly three zero modes
    u = np.repeat(np.sqrt(m), 3).reshape(-1, 3)                 # uniform translation, mass weighted
    for alpha in range(3):
        x = np.zeros_like(u)
        x[:, alpha] = u[:, alpha]
        assert np.abs(D @ x.ravel()).max() < 1e-13
    assert sla.eigsh(D, k=1, which="LA", return_eigenvectors=False)[0] > np.linalg.eigvalsh(
        hc.bloch(hc.POINTS["K"])).max()                         # the light defect has a mode above the band
