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


def bloch(k, p=DEFAULT, gauge="atom"):
    """6 x 6 mass-weighted dynamical matrix D(k), coordinates (A_x, A_y, A_z, B_x, B_y, B_z).

    gauge="atom": phases exp(ik.(r_b - r_a)) with the atomic positions; gauge="cell": phases with the
    cell vectors only, for which the torus matrix is a convolution over cells. Eigenvalues agree."""
    m = np.array([p["mA"], p["mB"]])
    D = np.zeros((6, 6), complex)
    for atoms, C in terms(p):
        for a, (sa, Ra) in enumerate(atoms):
            for b, (sb, Rb) in enumerate(atoms):
                d = position(sb, Rb) - position(sa, Ra)
                if gauge == "cell":
                    d = d - TAU[sb] + TAU[sa]
                D[3 * sa:3 * sa + 3, 3 * sb:3 * sb + 3] += C[a][b] * np.exp(1j * (k @ d)) / np.sqrt(m[sa] * m[sb])
    return D


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


def supercell(L, p=DEFAULT, defect=None):
    """Mass-weighted dynamical matrix of an L x L torus of primitive cells, 3 coordinates per atom.

    defect = (sublattice, (n1, n2), mass, factor) replaces one atom: its mass, and every energy term
    that contains it is multiplied by factor. Returns (D as CSR, masses per atom); atom index
    2 * (n1 * L + n2) + s, coordinate index 3 * atom + alpha."""
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
    return (S @ Phi @ S).tocsr(), m


def rigid_motions(points):
    """The six rigid motions (3 translations, 3 rotations) of atoms at the given in-plane positions."""
    r = np.c_[points, np.zeros(len(points))]
    out = [np.tile(e, (len(r), 1)) for e in np.eye(3)]
    out += [np.cross(axis, r) for axis in np.eye(3)]
    return out
