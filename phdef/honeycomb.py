"""Second toy model: a honeycomb lattice with two atoms per cell (like h-BN), displacements in x, y and z.

Bond length 1, lattice vectors a1 = (sqrt(3), 0) and a2 = (sqrt(3)/2, 3/2); atom A at the origin and
atom B at (0, 1). The harmonic energy is a sum of terms, each invariant under rigid translations and
rigid rotations:

* nearest-neighbour central springs k1 on the A-B bonds,
* second-neighbour central springs k2A (A-A) and k2B (B-B),
* an out-of-plane bending term for every atom, (kappa/2) (u_z,i - mean of u_z over its 3 neighbours)^2.

Hence the acoustic and rotational sum rules hold exactly, the flexural branch is quadratic
(omega_ZA ~ q^2), and in-plane and out-of-plane motion decouple (mirror plane z -> -z). A substitution
changes the mass of one atom and multiplies every term that contains that atom by a common factor.
"""
import numpy as np
import scipy.sparse as sp

from .chain import block_lanczos

S3 = np.sqrt(3.0)
A1, A2 = np.array([S3, 0.0]), np.array([S3 / 2, 1.5])
B1, B2 = 2 * np.pi * np.array([1 / S3, -1 / 3]), 2 * np.pi * np.array([0.0, 2 / 3])
TAU = (np.array([0.0, 0.0]), np.array([0.0, 1.0]))
NEIGHBOURS = {0: [(1, (0, 0)), (1, (1, -1)), (1, (0, -1))],   # the B neighbours of A in cell 0
              1: [(0, (0, 0)), (0, (-1, 1)), (0, (0, 1))]}    # the A neighbours of B in cell 0
SECOND = [(1, 0), (0, 1), (-1, 1)]                              # second-neighbour bonds, each once
POINTS = {"Γ": np.zeros(2), "M": B1 / 2, "K": (2 * B1 + B2) / 3}
DEFAULT = dict(mA=1.0, mB=1.3, k1=1.0, k2A=0.2, k2B=0.2, kapA=0.3, kapB=0.3)
IN_PLANE, OUT_OF_PLANE = [0, 1, 3, 4], [2, 5]                  # coordinate indices of the 6 x 6 blocks


def position(s, cell):
    """In-plane position of atom s (0 = A, 1 = B) in the cell (n1, n2)."""
    return TAU[s] + cell[0] * A1 + cell[1] * A2


def _central(k, a, b):
    e = np.r_[position(*b) - position(*a), 0.0]
    P = k * np.outer(e, e) / (e @ e)
    return [a, b], [[P, -P], [-P, P]]


def terms(p=DEFAULT):
    """Energy terms of one primitive cell: (atoms [(sublattice, cell)], blocks C[a][b] (3 x 3)) with
    E = (1/2) sum_ab u_a^T C[a][b] u_b."""
    out = [_central(p["k1"], (0, (0, 0)), nb) for nb in NEIGHBOURS[0]]
    out += [_central(p[k], (s, (0, 0)), (s, R)) for s, k in ((0, "k2A"), (1, "k2B")) for R in SECOND]
    zz = np.diag([0.0, 0.0, 1.0])
    c = np.array([1.0, -1 / 3, -1 / 3, -1 / 3])
    for s, kap in ((0, p["kapA"]), (1, p["kapB"])):
        out.append(([(s, (0, 0))] + NEIGHBOURS[s], [[kap * c[a] * c[b] * zz for b in range(4)] for a in range(4)]))
    return out


def bloch(k, p=DEFAULT, gauge="atom", polar=None):
    """6 x 6 mass-weighted dynamical matrix D(k), coordinates (A_x, A_y, A_z, B_x, B_y, B_z).

    gauge="atom": phases exp(ik.(r_b - r_a)) with the atomic positions; gauge="cell": phases with the
    cell vectors only, for which the torus matrix is a convolution over cells. Eigenvalues agree.
    polar: parameters of the long-range dipole-dipole term (see `long_range`), or None."""
    m = np.array([p["mA"], p["mB"]])
    D = np.zeros((6, 6), complex)
    for atoms, C in terms(p):
        for a, (sa, Ra) in enumerate(atoms):
            for b, (sb, Rb) in enumerate(atoms):
                d = position(sb, Rb) - position(sa, Ra)
                if gauge == "cell":
                    d = d - TAU[sb] + TAU[sa]
                D[3 * sa:3 * sa + 3, 3 * sb:3 * sb + 3] += C[a][b] * np.exp(1j * (k @ d)) / np.sqrt(m[sa] * m[sb])
    if polar is not None:
        D += long_range(k, p, polar, gauge)
    return D


# ---- polar switch: in-plane Born charges and the 2D Coulomb interaction ----
POLAR = dict(Z=0.48, reff=2.0, lam=0.7)   # Born charge +Z (A) / -Z (B) in plane, 2D screening length, dipole width 1/lam
AREA = abs(A1[0] * A2[1] - A1[1] * A2[0])
_G = np.array([n1 * B1 + n2 * B2 for n1 in range(-4, 5) for n2 in range(-4, 5)])


def _phi_long_range(k, polar):
    """Unweighted long-range force constants in the atom gauge: dipoles of Gaussian width 1/lam
    interacting through the screened 2D Coulomb kernel 2 pi / (|k| (1 + r_eff |k|)), summed over G."""
    Zs = (polar["Z"], -polar["Z"])
    kk = k[None, :] + _G
    kn = np.linalg.norm(kk, axis=1)
    w = np.zeros(len(kk))
    ok = kn > 1e-12                                   # the k -> 0 limit of k k / |k| vanishes in 2D
    w[ok] = np.exp(-kn[ok] ** 2 / (4 * polar["lam"] ** 2)) / (kn[ok] * (1 + polar["reff"] * kn[ok]))
    out = np.zeros((6, 6), complex)
    for s in range(2):
        for t in range(2):
            ph = w * np.exp(-1j * (_G @ (TAU[t] - TAU[s])))        # Poisson summation: phase e^{-iG.(tau_t - tau_s)}
            out[3 * s:3 * s + 2, 3 * t:3 * t + 2] = (2 * np.pi / AREA) * Zs[s] * Zs[t] * np.einsum("g,ga,gb->ab", ph, kk, kk)
    return out


def long_range(k, p=DEFAULT, polar=POLAR, gauge="atom"):
    """Mass-weighted long-range part of D(k). The q = 0 row sums are removed on site (acoustic sum
    rule); the in-plane Born charges leave out-of-plane motion untouched (mirror symmetry kept).
    In 2D the LO-TO splitting vanishes at Γ and grows linearly with |k|."""
    Phi = _phi_long_range(np.asarray(k, float), polar)
    P0 = _phi_long_range(np.zeros(2), polar)
    for s in range(2):
        Phi[3 * s:3 * s + 3, 3 * s:3 * s + 3] -= P0[3 * s:3 * s + 3, 0:3] + P0[3 * s:3 * s + 3, 3:6]
    if gauge == "cell":
        for s in range(2):
            for t in range(2):
                Phi[3 * s:3 * s + 3, 3 * t:3 * t + 3] *= np.exp(-1j * (k @ (TAU[t] - TAU[s])))
    m = np.repeat([p["mA"], p["mB"]], 3)
    return Phi / np.sqrt(np.outer(m, m))


def path(points=("Γ", "M", "K", "Γ"), per_unit=60):
    """k-points along a path through POINTS; returns (path coordinate, k (n, 2), tick positions)."""
    ks, ss, ticks = [], [], [0.0]
    for n, (a, b) in enumerate(zip(points[:-1], points[1:])):
        ka, kb = POINTS[a], POINTS[b]
        length = np.linalg.norm(kb - ka)
        t = np.linspace(0.0, 1.0, int(round(per_unit * length)) + 1)
        if n < len(points) - 2:
            t = t[:-1]
        ks.append(ka + np.outer(t, kb - ka))
        ss.append(ticks[-1] + length * t)
        ticks.append(ticks[-1] + length)
    return np.concatenate(ss), np.vstack(ks), ticks


def site_index(s, n1, n2, L):
    return 2 * ((n1 % L) * L + (n2 % L)) + s


def supercell(L, p=DEFAULT, defect=None, polar=None):
    """Mass-weighted dynamical matrix of an L x L torus of primitive cells, 3 coordinates per atom.

    defect = (sublattice, (n1, n2), mass, factor) replaces one atom: its mass, and every short-range
    energy term that contains it is multiplied by factor (its Born charge is kept). Returns (D,
    masses per atom); D is CSR, or a dense array with the long-range term when polar is given.
    Atom index 2 * (n1 * L + n2) + s, coordinate index 3 * atom + alpha."""
    if L < 3:
        raise ValueError("L must be at least 3")
    n1, n2 = (a.ravel() for a in np.meshgrid(np.arange(L), np.arange(L), indexing="ij"))
    d_site = site_index(defect[0], *defect[1], L) if defect else -1
    rows, cols, vals = [], [], []
    for atoms, C in terms(p):
        ids = [site_index(s, n1 + c[0], n2 + c[1], L) for s, c in atoms]
        factor = np.ones(L * L)
        if defect:
            factor[np.any(np.array(ids) == d_site, axis=0)] = defect[3]
        for a in range(len(atoms)):
            for b in range(len(atoms)):
                for al, be in zip(*np.nonzero(C[a][b])):
                    rows.append(3 * ids[a] + al)
                    cols.append(3 * ids[b] + be)
                    vals.append(factor * C[a][b][al, be])
    nat = 2 * L * L
    Phi = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                        shape=(3 * nat, 3 * nat)).tocsr()
    m = np.tile([p["mA"], p["mB"]], L * L)
    if defect:
        m[d_site] = defect[2]
    S = sp.diags(np.repeat(1.0 / np.sqrt(m), 3))
    if polar is None:
        return (S @ Phi @ S).tocsr(), m
    w = np.repeat(1.0 / np.sqrt(m), 3)                  # dense: the dipole interaction couples all atoms
    return (Phi.toarray() + torus_long_range(L, p, polar)) * np.outer(w, w), m


def torus_long_range(L, p=DEFAULT, polar=POLAR):
    """Unweighted long-range force constants of the L x L torus (dense, 6 L^2 x 6 L^2), the inverse
    Fourier transform of the cell-gauge long-range D(q) on the L x L mesh."""
    mvec = np.repeat([p["mA"], p["mB"]], 3)
    scale = np.sqrt(np.outer(mvec, mvec))
    m1, m2 = np.meshgrid(np.arange(L), np.arange(L), indexing="ij")
    Phq = np.array([long_range((a * B1 + b * B2) / L, p, polar, "cell") * scale for a, b in zip(m1.ravel(), m2.ravel())])
    PhR = np.fft.fft2(Phq.reshape(L, L, 6, 6), axes=(0, 1)).real / L ** 2     # Phi(s, t, R): sum_q e^{-iq.R}
    n1, n2 = (a.ravel() for a in np.meshgrid(np.arange(L), np.arange(L), indexing="ij"))
    blocks = PhR[(n1[None, :] - n1[:, None]) % L, (n2[None, :] - n2[:, None]) % L]  # [cell a, cell b] -> Phi(R_b - R_a)
    return blocks.transpose(0, 2, 1, 3).reshape(6 * L * L, 6 * L * L)


def defect_coupling(B0, coords, m_host, m_def):
    """Coupling block of the chain for a defect that changes masses inside the cluster: the columns
    of B0 (cluster coordinates) scale with sqrt(m_host / m_def). Exact as long as the force constants
    between cluster and bath are unchanged; needed when the bath couples to the defect atom directly,
    as with long-range forces."""
    atom = np.asarray(coords) // 3
    return B0 * np.sqrt(m_host[atom] / m_def[atom])[None, :]


def atom_neighbours(atom, L):
    """Indices of the three nearest neighbours of an atom of the L x L torus."""
    c, s = divmod(atom, 2)
    n1, n2 = divmod(c, L)
    return [site_index(t, n1 + R[0], n2 + R[1], L) for t, R in NEIGHBOURS[s]]


def hop_distance(site, L):
    """Number of nearest-neighbour bonds from atom `site` to every atom of the L x L torus."""
    dist = np.full(2 * L * L, -1)
    dist[site] = 0
    front = [site]
    while front:
        nxt = []
        for a in front:
            for b in atom_neighbours(a, L):
                if dist[b] < 0:
                    dist[b] = dist[a] + 1
                    nxt.append(b)
        front = nxt
    return dist


def coordinates(atoms):
    """Coordinate indices (3 per atom) of a list of atoms."""
    return (3 * np.asarray(atoms)[:, None] + np.arange(3)[None, :]).ravel()


def bloch_modes(L, p=DEFAULT, polar=None):
    """Host modes on the L x L q-mesh: omega^2 (6 L^2,) and eigenvectors eps[q, coordinate, branch]
    of the cell-gauge Bloch matrix, q = (m1 b1 + m2 b2) / L with q index m1 * L + m2."""
    m1, m2 = (a.ravel() for a in np.meshgrid(np.arange(L), np.arange(L), indexing="ij"))
    qs = (np.outer(m1, B1) + np.outer(m2, B2)) / L
    w2, eps = np.linalg.eigh(np.array([bloch(q, p, gauge="cell", polar=polar) for q in qs]))
    return w2, eps, np.c_[m1, m2]


def mode_start_block(coords, L, eps, mq):
    """<q nu | e_c> for the real-space coordinates e_c: rows are the modes (q-major, then branch)."""
    atom, alpha = np.divmod(np.asarray(coords), 3)
    cell, s = np.divmod(atom, 2)
    n1, n2 = np.divmod(cell, L)
    phase = np.exp(-2j * np.pi * (np.outer(mq[:, 0], n1) + np.outer(mq[:, 1], n2)) / L)   # e^{-iq.R_c}
    comp = eps[:, 3 * s + alpha, :].conj()                                                # (nq, nc, 6)
    return (phase[:, :, None] * comp / L).transpose(0, 2, 1).reshape(-1, len(coords))


def mode_chain(coords, L, m, p=DEFAULT, polar=None):
    """Chain of the L x L host bath built from the Bloch modes (the reciprocal-space route).

    Block Lanczos on diag(omega^2) of the 6 L^2 host modes, started from the cluster coordinates
    `coords`: the first block is the cluster itself, the next ones are the chain of the bath, so no
    real-space bath matrix is needed. Returns (B0, A, B) for `upfold.assemble` with B0 acting on the
    cluster coordinates, and (A1, W): the first block A1 = W^H D0_cc W with the cluster rotation W."""
    w2, eps, mq = bloch_modes(L, p, polar)
    S = mode_start_block(coords, L, eps, mq)
    B0f, Af, Bf = block_lanczos(lambda X: w2.ravel()[:, None] * X, S, m + 1)
    W = B0f.conj().T                        # S = Q1 B0f with S orthonormal, so Q1 = S W
    return Bf[0] @ W.conj().T, Af[1:], Bf[1:], (Af[0], W)


def rigid_motions(points):
    """The six rigid motions (3 translations, 3 rotations) of atoms at the given in-plane positions."""
    r = np.c_[points, np.zeros(len(points))]
    out = [np.tile(e, (len(r), 1)) for e in np.eye(3)]
    out += [np.cross(axis, r) for axis in np.eye(3)]
    return out


# ---- one substitution in displacement coordinates: T-matrix and averaged self-energy ----
def masses(p=DEFAULT):
    """Masses of the six coordinates (A_x, A_y, A_z, B_x, B_y, B_z)."""
    return np.repeat([p["mA"], p["mB"]], 3)


def perturbation(defect, Lb=8, p=DEFAULT, polar=None, radius=2):
    """Support of a substitution (atoms within `radius` bonds of the A defect in cell 0): cell offsets,
    sublattices, and the change of the unweighted force constants dPhi (3n x 3n) and of the masses dM
    (3n), taken from the real-space supercell; asserts that nothing changes outside the support."""
    site = site_index(0, 0, 0, Lb)
    atoms = np.where(hop_distance(site, Lb) <= radius)[0]
    cell, s = np.divmod(atoms, 2)
    n1, n2 = np.divmod(cell, Lb)
    offs = np.c_[np.where(n1 > Lb // 2, n1 - Lb, n1), np.where(n2 > Lb // 2, n2 - Lb, n2)]
    D0, m0 = supercell(Lb, p, polar=polar)
    D, m = supercell(Lb, p, defect=defect, polar=polar)
    D0, D = (X.toarray() if sp.issparse(X) else X for X in (D0, D))
    Phi0 = D0 * np.outer(np.repeat(np.sqrt(m0), 3), np.repeat(np.sqrt(m0), 3))
    Phi = D * np.outer(np.repeat(np.sqrt(m), 3), np.repeat(np.sqrt(m), 3))
    C = coordinates(atoms)
    rest = np.setdiff1d(np.arange(len(D)), C)
    assert np.abs((Phi - Phi0)[rest]).max() < 1e-12
    return offs, s, (Phi - Phi0)[np.ix_(C, C)], np.repeat(m[atoms] - m0[atoms], 3)


def host_green_support(omegas, eta, Lq, offs, subl, p=DEFAULT, polar=None):
    """Displacement Green function of the host between the support coordinates, (nw, 3n, 3n), from the Bloch
    modes of an Lq x Lq mesh: g0(R_a - R_b) = (1/N) sum_q e^{iq.(R_a - R_b)} [z M - Phi(q)]^-1, by FFT."""
    w2, eps, _ = bloch_modes(Lq, p, polar)
    U = eps / np.sqrt(masses(p))[None, :, None]
    n = len(offs)
    d = (offs[:, None, :] - offs[None, :, :]) % Lq
    idx = 3 * np.asarray(subl)[:, None] + np.arange(3)[None, :]
    out = np.empty((len(omegas), 3 * n, 3 * n), complex)
    for i, w in enumerate(omegas):
        F = np.einsum("qik,qk,qjk->qij", U, 1.0 / ((w + 1j * eta) ** 2 - w2), U.conj()).reshape(Lq, Lq, 6, 6)
        g = np.fft.ifft2(F, axes=(0, 1))[d[:, :, 0], d[:, :, 1]]            # (n, n, 6, 6)
        out[i] = g[np.arange(n)[:, None, None, None], np.arange(n)[None, :, None, None],
                   idx[:, None, :, None], idx[None, :, None, :]].transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    return out


def tmatrix_u(omegas, eta, g0, dPhi, dM):
    """T = (1 - V g0)^-1 V with V(z) = dPhi - z dM (displacement coordinates)."""
    z = (np.asarray(omegas) + 1j * eta) ** 2
    V = dPhi[None] - z[:, None, None] * np.diag(dM)[None]
    return np.linalg.solve(np.eye(len(dM))[None] - V @ g0, V)


def self_energy(q, T, offs, subl, nd):
    """Sigma_st(q) = n_d sum_{a in s, b in t} e^{-iq.(R_a - R_b)} T_ab, (nw, 6, 6), in the convention of
    bloch(gauge="cell"), D_st(q) = sum_R D(0s, Rt) e^{iq.R}. The lattice has no inversion centre, so the
    sign of the phase matters: Sigma(-q) = Sigma(q)^T."""
    ph = np.exp(1j * (offs @ np.array([q @ A1, q @ A2])))                   # <R_a s|q s> = e^{iq.R_a}
    P = np.zeros((3 * len(offs), 6), complex)
    for a in range(len(offs)):
        P[3 * a:3 * a + 3, 3 * subl[a]:3 * subl[a] + 3] = ph[a] * np.eye(3)
    return nd * np.einsum("ai,wab,bj->wij", P.conj(), T, P)
