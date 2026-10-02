"""Block Lanczos chain of a Hermitian operator, started from a coupling block.

Used to upfold the host phonon bath seen by a defect cluster: run on the bulk
dynamical matrix D_BB, started from the cluster->bulk coupling D_BC = Q1 @ B0.
The chain blocks (A_j, B_j) are independent of the frequency and -- as long as
the cluster contains every mass and force-constant change -- of the defect.
"""
import numpy as np


def _herm(X):
    return X.conj().T


def block_lanczos(op, start, m, tol=1e-10, reorth="full", return_Q=False):
    """Block-tridiagonalise a Hermitian operator on the block Krylov space of ``start``.

    Parameters
    ----------
    op : (n, n) ndarray, scipy sparse matrix, or callable ``X -> op @ X``
        Hermitian operator (here: the host bulk dynamical matrix D_BB).
    start : (n, p) ndarray
        Source block (here: the coupling D_BC from the cluster into the bulk).
    m : int
        Maximum number of Lanczos blocks.
    tol : float
        Relative singular-value cutoff used for rank revealing and deflation.
        The recursion stops early if the Krylov space becomes invariant.
    reorth : {"full", "local"}
        "full" re-orthogonalises each new block against all previous blocks (twice);
        robust, O(n (m r)^2) work -- fine for prototypes. "local" keeps only the
        three-term recurrence (cheap; loses orthogonality on long chains).
    return_Q : bool
        Also return the list of Lanczos blocks Q_j (needed e.g. to project vertices
        onto chain sites).

    Returns
    -------
    B0 : (r0, p) ndarray
        ``start = Q1 @ B0``.
    A : list of (r_j, r_j) ndarrays
        Diagonal blocks ``Q_j^H op Q_j``.
    B : list of (r_{j+1}, r_j) ndarrays
        Couplings ``Q_{j+1}^H op Q_j``; ``len(B) == len(A) - 1``.
    Q : list of (n, r_j) ndarrays, only if ``return_Q``.
    """
    apply = op if callable(op) else (lambda X: op @ X)
    if reorth not in ("full", "local"):
        raise ValueError("reorth must be 'full' or 'local'")

    U, s, Vh = np.linalg.svd(np.asarray(start), full_matrices=False)
    if s.size == 0 or s[0] == 0.0:
        raise ValueError("start block is zero")
    r = int((s > tol * s[0]).sum())
    Qj = U[:, :r]
    B0 = s[:r, None] * Vh[:r]

    A, B, Q = [], [], [Qj]
    Qprev, Bprev = None, None
    scale = s[0]
    for _ in range(m):
        W = np.asarray(apply(Qj))
        Aj = _herm(Qj) @ W
        Aj = 0.5 * (Aj + _herm(Aj))
        A.append(Aj)
        if len(A) == m:
            break
        W = W - Qj @ Aj
        if Qprev is not None:
            W = W - Qprev @ _herm(Bprev)
        if reorth == "full":
            Qall = np.hstack(Q)
            for _pass in range(2):
                W = W - Qall @ (_herm(Qall) @ W)
        scale = max(scale, np.linalg.norm(Aj, 2))
        U, s, Vh = np.linalg.svd(W, full_matrices=False)
        r = int((s > tol * scale).sum())
        if r == 0:                      # invariant Krylov space reached: chain is exact
            break
        Qn = U[:, :r]
        Bj = s[:r, None] * Vh[:r]       # W = Qn @ Bj
        B.append(Bj)
        Q.append(Qn)
        Qprev, Bprev, Qj = Qj, Bj, Qn
    if return_Q:
        return B0, A, B, Q[:len(A)]
    return B0, A, B
