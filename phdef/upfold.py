"""Static upfolded dynamical matrix: defect cluster + Krylov chain of the host bath.

Lihm & Park (arXiv:2409.07393, SM S14-S18) turn a frequency-dependent phonon
self-energy into a frequency-independent matrix by adding auxiliary modes.
For a point defect the frequency dependence comes from folding the infinite
host onto the defect cluster C,

    D_CC(z) = [z - D_CC - Pi_B(z)]^-1,   Pi_B(z) = D_CB (z - D_BB)^-1 D_BC,
    z = (omega + i eta)^2,

and the auxiliary modes are the sites of a block Lanczos chain of D_BB:

    Ct = [[D_CC,   B0^H E1^H],
          [E1 B0,  L_m      ]],      D_CC(z) = [(z - Ct)^-1]_CC .

Everything here is frequency independent; one diagonalisation of Ct gives the
cluster Green function at every frequency (sum over "dressed" modes).
"""
import numpy as np


def _herm(X):
    return X.conj().T


def assemble(D_CC, B0, A, B):
    """Assemble Ct from the cluster block and a chain from :func:`block_lanczos`."""
    D_CC = np.asarray(D_CC)
    nC = D_CC.shape[0]
    sizes = [a.shape[0] for a in A]
    off = np.concatenate([[nC], nC + np.cumsum(sizes)]).astype(int)
    dtype = np.result_type(D_CC, B0, *A, *B) if B else np.result_type(D_CC, B0, *A)
    Ct = np.zeros((off[-1], off[-1]), dtype=dtype)
    Ct[:nC, :nC] = D_CC
    r0 = sizes[0]
    Ct[nC:nC + r0, :nC] = B0
    Ct[:nC, nC:nC + r0] = _herm(B0)
    for j, Aj in enumerate(A):
        o, r = off[j], sizes[j]
        Ct[o:o + r, o:o + r] = Aj
        if j < len(B):
            o2, r2 = off[j + 1], sizes[j + 1]
            Ct[o2:o2 + r2, o:o + r] = B[j]
            Ct[o:o + r, o2:o2 + r2] = _herm(B[j])
    return Ct


def modes(Ct):
    """Dressed modes of the upfolded matrix: eigenvalues omega^2 and eigenvectors."""
    return np.linalg.eigh(Ct)


def cluster_green_modes(w2, V, nC, z):
    """G_CC(z) = sum_beta e_C,beta e_C,beta^H / (z - omega_beta^2) for an array of z.

    Returns an array of shape (len(z), nC, nC).
    """
    z = np.atleast_1d(z)
    VC = V[:nC]
    return np.einsum("ib,zb,jb->zij", VC, 1.0 / (z[:, None] - w2[None, :]), VC.conj())


def cluster_green_cf(D_CC, B0, A, B, z):
    """Same G_CC(z) by the block continued fraction (no diagonalisation; cost per z).

    Evaluated bottom-up: S_m = (z - A_m)^-1, S_j = (z - A_j - B_j^H S_{j+1} B_j)^-1,
    G_CC = (z - D_CC - B0^H S_1 B0)^-1.
    """
    z = np.atleast_1d(z)
    nC = np.asarray(D_CC).shape[0]
    out = np.empty((z.size, nC, nC), dtype=complex)
    for iz, zz in enumerate(z):
        S = np.linalg.inv(zz * np.eye(A[-1].shape[0]) - A[-1])
        for j in range(len(A) - 2, -1, -1):
            S = np.linalg.inv(zz * np.eye(A[j].shape[0]) - A[j] - _herm(B[j]) @ S @ B[j])
        out[iz] = np.linalg.inv(zz * np.eye(nC) - D_CC - _herm(B0) @ S @ B0)
    return out


def damped_companion(Ct, gamma):
    """Lihm-Park Eq. 9 / S23: K = [[-i G, 1], [Ct, -i G]] with G = diag(gamma).

    [(omega + i G)^2 - Ct]^-1 = (1 0) (omega - K)^-1 (0 1)^T, and the complex
    eigenvalues of K, omega~ - i gamma~, are the damped mode frequencies.
    """
    n = Ct.shape[0]
    G = np.diag(np.asarray(gamma, dtype=float) * np.ones(n))
    return np.block([[-1j * G, np.eye(n)], [Ct, -1j * G]])


def green_from_companion(K, omegas):
    """[(omega + i G)^2 - Ct]^-1 from one non-Hermitian eigendecomposition of K."""
    n = K.shape[0] // 2
    lam, R = np.linalg.eig(K)
    L = np.linalg.inv(R)
    omegas = np.atleast_1d(omegas)
    return np.einsum("ib,wb,bj->wij", R[:n], 1.0 / (omegas[:, None] - lam[None, :]), L[:, n:])
