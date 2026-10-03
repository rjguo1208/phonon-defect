"""Tests of the honeycomb toy model (run: python -m pytest tests)."""
import os
import sys

import numpy as np
import scipy.sparse.linalg as sla

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
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


def _cluster(L, radius=2):
    site = hc.site_index(0, 0, 0, L)
    hop = hc.hop_distance(site, L)
    return site, hc.coordinates(np.where(hop <= radius)[0]), hc.coordinates(np.where(hop > radius)[0])


def test_cluster_of_two_bonds_contains_every_changed_term():
    L = 8
    site, C, B = _cluster(L)
    assert len(C) == 30                                          # 1 + 3 + 6 atoms
    D0, _ = hc.supercell(L)
    for d in [(0, (0, 0), 0.3, 1.3), (0, (0, 0), 8.0, 0.6)]:
        D, _ = hc.supercell(L, defect=d)
        assert abs(D[B][:, C] - D0[B][:, C]).max() < 1e-14      # coupling block is a host quantity
        assert abs(D[B][:, B] - D0[B][:, B]).max() < 1e-14      # and so is the bath


def test_mode_chain_reproduces_the_defect_supercell():
    """On the same q-mesh the chain terminates, so the upfolded matrix equals the supercell."""
    from phdef import assemble, modes, cluster_green_modes
    L = 8
    site, C, B = _cluster(L)
    B0, A, Bc, (A1, W) = hc.mode_chain(C, L, 500)
    D0, _ = hc.supercell(L)
    assert np.abs(A1 - W.conj().T @ D0[C][:, C].toarray() @ W).max() < 1e-12
    z = (np.array([0.2, 0.9, 1.45, 2.0]) + 0.03j) ** 2
    for d in [None, (0, (0, 0), 0.3, 1.3), (0, (0, 0), 8.0, 0.6)]:
        D = D0 if d is None else hc.supercell(L, defect=d)[0]
        G = cluster_green_modes(*modes(assemble(D[C][:, C].toarray(), B0, A, Bc)), len(C), z)
        Dd = D.toarray()
        Gx = np.array([np.linalg.inv(zz * np.eye(len(Dd)) - Dd)[np.ix_(C, C)] for zz in z])
        assert np.abs(G - Gx).max() < 1e-9 * np.abs(Gx).max()


def test_polar_switch_sum_rules_and_lo_to():
    P = hc.POLAR
    w2 = np.linalg.eigvalsh(hc.bloch(hc.POINTS["Γ"], polar=P))
    assert np.abs(w2[:3]).max() < 1e-12                          # three zero modes at Γ
    assert abs(w2[-1] - w2[-2]) < 1e-12                           # in 2D, LO = TO at Γ
    split = [np.diff(np.sqrt(np.linalg.eigvalsh(hc.bloch(np.array([q, 0.0]), polar=P)))[-2:])[0] / q
             for q in (1e-3, 2e-3)]
    assert split[0] > 0.2 and abs(split[0] - split[1]) < 1e-2 * split[0]   # the splitting grows linearly
    D = hc.bloch(np.array([0.3, 0.1]), polar=P)
    assert np.abs(D[np.ix_(hc.IN_PLANE, hc.OUT_OF_PLANE)]).max() < 1e-14     # mirror symmetry kept
    n = 24
    assert min(np.linalg.eigvalsh(hc.bloch((i * hc.B1 + j * hc.B2) / n, polar=P)).min()
               for i in range(n) for j in range(n)) > -1e-12


def test_polar_torus_equals_bloch_spectrum():
    L = 5
    D, _ = hc.supercell(L, polar=hc.POLAR)
    grid = np.sort(np.concatenate([np.linalg.eigvalsh(hc.bloch((i * hc.B1 + j * hc.B2) / L, polar=hc.POLAR))
                                   for i in range(L) for j in range(L)]))
    assert np.allclose(np.linalg.eigvalsh(D), grid, atol=1e-11)


def test_polar_mass_defect_needs_the_column_scaling():
    """With long-range forces the bath couples to the defect atom itself, so the coupling block
    changes with its mass; scaling the columns of B0 restores the exact result."""
    from phdef import assemble, modes, cluster_green_modes
    L = 8
    site, C, B = _cluster(L)
    B0, A, Bc, _ = hc.mode_chain(C, L, 500, polar=hc.POLAR)
    D0, m0 = hc.supercell(L, polar=hc.POLAR)
    D, m = hc.supercell(L, defect=(0, (0, 0), 0.3, 1.3), polar=hc.POLAR)
    z = (np.array([0.5, 1.45, 2.0]) + 0.03j) ** 2
    Gx = np.array([np.linalg.inv(zz * np.eye(len(D)) - D)[np.ix_(C, C)] for zz in z])
    Dcc = D[np.ix_(C, C)]
    G = cluster_green_modes(*modes(assemble(Dcc, hc.defect_coupling(B0, C, m0, m), A, Bc)), len(C), z)
    G_naive = cluster_green_modes(*modes(assemble(Dcc, B0, A, Bc)), len(C), z)
    assert np.abs(G - Gx).max() < 1e-9 * np.abs(Gx).max()
    assert np.abs(G_naive - Gx).max() > 1e-4 * np.abs(Gx).max()


def test_position_average_of_one_defect_is_the_first_order_self_energy():
    """One defect on an L x L torus, averaged over its N positions: <g>(q) = g0 + g0 Sigma g0 / N exactly,
    with Sigma from the T-matrix on the support. Fixes the sign of the Bloch phase in Sigma(q)."""
    import honeycomb_figures as hf
    L, eta = 8, 0.05
    om = np.array([0.3, 1.0, 1.45, 2.0, 3.1])
    z = (om + 1j * eta) ** 2
    BM = np.c_[hc.B1, hc.B2]
    n1, n2 = (a.ravel() for a in np.meshgrid(np.arange(L), np.arange(L), indexing="ij"))
    R = np.outer(n1, hc.A1) + np.outer(n2, hc.A2)
    for d in [(0, (0, 0), 0.3, 1.3), (0, (0, 0), 8.0, 0.6)]:
        offs, subl, dPhi, dM = hf.perturbation(d)
        T = hf.tmatrix_u(om, eta, hf.host_green_support(om, eta, L, offs, subl), dPhi, dM)
        D, m = hc.supercell(L, defect=d)
        w = np.repeat(1.0 / np.sqrt(m), 3)
        G = np.array([w[:, None] * np.linalg.inv(zz * np.eye(len(w)) - D.toarray()) * w[None, :] for zz in z])
        for mq in [(1, 0), (2, 1), (5, 2), (1, 3)]:                  # none of them equivalent to -q
            q = BM @ np.array(mq, float) / L
            F = np.kron(np.exp(1j * (R @ q))[:, None], np.eye(6))   # columns: Bloch waves at q
            avg = np.einsum("ai,wab,bj->wij", F.conj(), G, F) / L ** 2
            g0 = np.linalg.inv(z[:, None, None] * np.diag(hf.MVEC)[None]
                               - (hc.bloch(q, gauge="cell") * np.sqrt(np.outer(hf.MVEC, hf.MVEC)))[None])
            first = g0 + g0 @ hf.self_energy(q, T, offs, subl, 1.0 / L ** 2) @ g0
            assert np.abs(first - avg).max() < 1e-10 * np.abs(avg).max()
