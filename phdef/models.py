"""Toy lattice models used by the examples and tests."""
import numpy as np
import scipy.sparse as sp


def square_lattice_dynmat(L, Mdef=1.0, Kdef=1.0):
    """Mass-weighted dynamical matrix of an L x L square lattice (scalar displacements).

    Nearest-neighbour springs K = 1 and masses M = 1 (open boundaries), with a point
    defect at the centre site: mass ``Mdef`` and spring ``Kdef`` on its four bonds.
    The on-site terms are built from the springs, so the acoustic sum rule holds.
    Weighting uses the *defect* masses, so the matrix is static and Hermitian.

    Returns (D, centre_index) with D a CSR matrix of size L*L.
    """
    N = L * L
    c0 = (L // 2) * L + L // 2
    rows, cols, vals = [], [], []
    diag = np.zeros(N)
    for i in range(L):
        for j in range(L):
            a = i * L + j
            for di, dj in ((1, 0), (0, 1)):
                ii, jj = i + di, j + dj
                if ii < L and jj < L:
                    b = ii * L + jj
                    k = Kdef if (a == c0 or b == c0) else 1.0
                    rows += [a, b]
                    cols += [b, a]
                    vals += [-k, -k]
                    diag[a] += k
                    diag[b] += k
    idx = np.arange(N)
    Phi = sp.coo_matrix((np.r_[vals, diag], (np.r_[rows, idx], np.r_[cols, idx])),
                        shape=(N, N)).tocsr()
    mass = np.ones(N)
    mass[c0] = Mdef
    s = sp.diags(1.0 / np.sqrt(mass))
    return (s @ Phi @ s).tocsr(), c0


def manhattan_partition(L, RC):
    """Cluster = sites within Manhattan distance RC of the centre; bulk = the rest."""
    ic = L // 2
    i, j = np.divmod(np.arange(L * L), L)
    dist = np.abs(i - ic) + np.abs(j - ic)
    return np.where(dist <= RC)[0], np.where(dist > RC)[0]
