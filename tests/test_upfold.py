"""Tests for the cluster + Krylov-chain upfolding (run: python -m pytest tests)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from phdef import (block_lanczos, assemble, modes, cluster_green_modes, cluster_green_cf,  # noqa: E402
                   damped_companion, green_from_companion)
from phdef.models import square_lattice_dynmat, manhattan_partition  # noqa: E402

L, RC = 41, 2


@pytest.fixture(scope="module")
def host_chain():
    C, B = manhattan_partition(L, RC)
    host, c0 = square_lattice_dynmat(L)
    D_BB = host[B][:, B].tocsr()
    D_BC = host[B][:, C].toarray()
    B0, A, Bc, Q = block_lanczos(D_BB, D_BC, 60, return_Q=True)
    return dict(C=C, B=B, c0=c0, D_BC=D_BC, B0=B0, A=A, Bc=Bc, Q=Q)


def exact_cluster_green(D, C, z):
    Dd = D.toarray()
    out = []
    for zz in z:
        G = np.linalg.inv(zz * np.eye(Dd.shape[0]) - Dd)
        out.append(G[np.ix_(C, C)])
    return np.array(out)


def test_lanczos_blocks_orthonormal(host_chain):
    Q = np.hstack(host_chain["Q"])
    assert np.allclose(Q.T @ Q, np.eye(Q.shape[1]), atol=1e-10)


def test_coupling_block_is_defect_independent(host_chain):
    C, B = host_chain["C"], host_chain["B"]
    for Md, Kd in [(0.3, 1.3), (8.0, 0.6)]:
        D, _ = square_lattice_dynmat(L, Md, Kd)
        assert np.allclose(D[B][:, C].toarray(), host_chain["D_BC"])


def test_local_mode_exact_with_short_chain(host_chain):
    """A local mode above the band is a bound state: a short chain already gives it exactly."""
    D, _ = square_lattice_dynmat(L, 0.3, 1.3)
    C = host_chain["C"]
    m = 12
    Ct = assemble(D[C][:, C].toarray(), host_chain["B0"], host_chain["A"][:m], host_chain["Bc"][:m - 1])
    w2, _ = modes(Ct)
    exact_top = np.linalg.eigvalsh(D.toarray())[-1]
    assert exact_top > 8.0                      # above the host band (omega_max^2 = 8)
    assert abs(w2[-1] - exact_top) < 1e-9


def test_continued_fraction_equals_mode_sum(host_chain):
    D, _ = square_lattice_dynmat(L, 8.0, 0.6)
    C = host_chain["C"]
    m = 30
    D_CC = D[C][:, C].toarray()
    A, Bc = host_chain["A"][:m], host_chain["Bc"][:m - 1]
    z = (np.array([0.3, 0.5, 1.0, 2.0]) + 0.05j) ** 2
    w2, V = modes(assemble(D_CC, host_chain["B0"], A, Bc))
    G_modes = cluster_green_modes(w2, V, len(C), z)
    G_cf = cluster_green_cf(D_CC, host_chain["B0"], A, Bc, z)
    assert np.max(np.abs(G_modes - G_cf)) < 1e-10 * np.max(np.abs(G_cf))


@pytest.mark.parametrize("Md,Kd", [(1.0, 1.0), (0.3, 1.3), (8.0, 0.6)])
def test_converges_to_exact_cluster_green(host_chain, Md, Kd):
    """Pristine (null test), light and heavy defect: same chain, only D_CC changes."""
    D, _ = square_lattice_dynmat(L, Md, Kd)
    C = host_chain["C"]
    D_CC = D[C][:, C].toarray()
    z = (np.array([0.3, 0.5, 1.0, 2.0, 2.7]) + 0.3j) ** 2
    G_exact = exact_cluster_green(D, C, z)
    G_cf = cluster_green_cf(D_CC, host_chain["B0"], host_chain["A"], host_chain["Bc"], z)
    assert np.max(np.abs(G_cf - G_exact)) < 1e-6 * np.max(np.abs(G_exact))


def test_damped_companion_matches_quadratic_inverse():
    rng = np.random.default_rng(3)
    n = 6
    X = rng.normal(size=(n, n))
    Ct = X @ X.T + n * np.eye(n)               # Hermitian, positive: squared frequencies
    gamma = rng.uniform(0.05, 0.5, size=n)
    K = damped_companion(Ct, gamma)
    omegas = np.array([0.4, 1.3, 2.2, 3.1])
    G = green_from_companion(K, omegas)
    for w, Gw in zip(omegas, G):
        a = w * np.eye(n) + 1j * np.diag(gamma)
        assert np.allclose(Gw, np.linalg.inv(a @ a - Ct), atol=1e-10)
    lam = np.linalg.eigvals(K)
    assert np.all(lam.imag < 0)                 # retarded: damped modes in the lower half plane
