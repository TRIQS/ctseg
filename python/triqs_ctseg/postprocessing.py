"""
Post-processing routines for the CT-SEG solver.

This module contains self-energy extraction (Dyson or improved estimator),
polarizability and screened interaction calculations, and related helpers.
"""

import triqs.utility.mpi as mpi
import numpy as np
from itertools import product

from triqs.gfs import (
    MeshImFreq, Gf, BlockGf,
    make_gf_dlr, make_gf_imfreq, make_hermitian
)
from triqs.gfs.tools import make_zero_tail
from triqs.operators.util.extractors import extract_U_dict2, dict_to_matrix

from triqs.solver_utils import tail_fit


# =============================================================================
# Helpers: color <-> block index mapping
# =============================================================================

def build_color_tables(gf_struct):
    """Build per-color lookup tables for a gf_struct.

    Returns (block_name, index_in_block, orbital), lists of length n_color.
    The orbital index counts up/down orbitals separately, so that each
    spin sector yields 0..n_orb-1.
    """
    block_name, index_in_block, orbital = [], [], []
    offsets = {"up": 0, "dn": 0}
    for blk_name, blk_dim in gf_struct:
        key = "up" if blk_name[:2] == "up" else "dn"
        for i in range(blk_dim):
            block_name.append(blk_name)
            index_in_block.append(i)
            orbital.append(offsets[key] + i)
        offsets[key] += blk_dim
    return block_name, index_in_block, orbital


# =============================================================================
# Interaction tensor extraction
# =============================================================================

def _dd_to_4idx(direct, exchange):
    """Expand a density-density matrix into a 4-index tensor.

    direct[i, j]   -> entries of type (i, j, i, j)
    exchange[i, j] -> entries of type (i, j, j, i) and (i, i, j, j), for i != j
    """
    n = direct.shape[0]
    out = np.zeros((n, n, n, n), dtype=complex)
    for i, j in product(range(n), repeat=2):
        out[i, j, i, j] = direct[i, j]
        if i != j:
            out[i, j, j, i] = exchange[i, j]
            out[i, i, j, j] = exchange[i, j]
    return out


def _dd_matrix_to_4idx(M):
    """Split a 2n_orb x 2n_orb density-density matrix (up/down layout) into
    direct (opposite-spin) and exchange (direct - same-spin) blocks, then
    expand to a 4-index tensor."""
    n_orb = M.shape[0] // 2
    direct = M[:n_orb, n_orb:2*n_orb]
    exchange = direct - M[:n_orb, :n_orb]
    return _dd_to_4idx(direct, exchange)


def extract_u_tensor_from_h_int(h_int, gf_struct, return_4idx=False):
    """
    Return U tensor in density-density basis from a Coulomb many-body operator h_int.
    If return_4idx is True, construct the full 4-index Uijkl tensor.
    """
    U_dd = dict_to_matrix(extract_U_dict2(h_int), gf_struct=gf_struct)
    return _dd_matrix_to_4idx(U_dd) if return_4idx else U_dd


# =============================================================================
# Spectral truncation
# =============================================================================

def check_spectrum(A, radius=1.0, truncation=False):
    """Check eigenvalues of A and optionally truncate those with |lambda| > radius."""
    A = np.asarray(A)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("A must be square.")

    w, V = np.linalg.eig(A)
    keep = np.abs(w) < radius
    w_drop = w[~keep]

    if w_drop.size > 0:
        mpi.report(f"Unqualified eigenvalues (|lambda| > {radius}):")
        for lam in w_drop:
            mpi.report(f"  {lam}")
        mpi.report("")

    if not truncation or w_drop.size == 0:
        return A

    mpi.report(f"Reconstructing matrix with eigenvalues |lambda| <= {radius} only.\n")
    w_keep = w[keep]
    V_keep = V[:, keep]
    if w_keep.size == 0:
        return np.zeros_like(A, dtype=np.result_type(A, 1j))
    return V_keep @ np.diag(w_keep) @ np.linalg.pinv(V_keep)


# =============================================================================
# Screening matrix extraction
# =============================================================================

def _assemble_color_gf(src_block2gf, block_name, index_in_block):
    """Assemble a (n_color, n_color)-target Gf from a Block2Gf by picking
    the (index_in_block[c1], index_in_block[c2]) slot of block
    (block_name[c1], block_name[c2]) for each color pair."""
    n_color = len(block_name)
    mesh = src_block2gf[block_name[0], block_name[0]].mesh
    out = Gf(mesh=mesh, target_shape=(n_color, n_color))
    for c1, c2 in product(range(n_color), repeat=2):
        out.data[:, c1, c2] = (
            src_block2gf[block_name[c1], block_name[c2]].data[:, index_in_block[c1], index_in_block[c2]]
        )
    return out


def extract_screen_matrix_from_D0_tau(blk2_D0_tau, gf_struct, return_4idx=False):
    block_name, index_in_block, _ = build_color_tables(gf_struct)
    D0_tau = _assemble_color_gf(blk2_D0_tau, block_name, index_in_block)

    w0_mesh = MeshImFreq(beta=D0_tau.mesh.beta, statistic="Boson", n_iw=1)
    D0_iw = Gf(mesh=w0_mesh, target_shape=D0_tau.target_shape)
    D0_iw.set_from_fourier(D0_tau, make_zero_tail(D0_iw, n_moments=2))
    Dw0_dd = D0_iw.data[0].real

    return _dd_matrix_to_4idx(Dw0_dd) if return_4idx else Dw0_dd


# =============================================================================
# Hartree-Fock self-energy
# =============================================================================

def compute_sigma_hartreefock(solver):
    """
    Compute the static Hartree-Fock self-energy with density-density interactions.
    """
    mpi.report('\nEvaluating static impurity self-energy analytically in density-density basis:')

    block_name, index_in_block, _ = build_color_tables(solver.gf_struct)
    n_color = len(block_name)

    V = extract_u_tensor_from_h_int(h_int=solver.h_int, gf_struct=solver.gf_struct)
    if solver.D0_tau is not None:
        V += extract_screen_matrix_from_D0_tau(blk2_D0_tau=solver.D0_tau, gf_struct=solver.gf_struct)

    densities = np.array([
        solver.results.densities[block_name[c]][index_in_block[c]]
        for c in range(n_color)
    ], dtype=float)

    # For each color c1, Sigma_HF_diag[c1] = sum_c2 V[c1,c2] * densities[c2]
    hf_per_color = V.real @ densities

    Sigma_HartreeFock = {}
    c = 0
    for blk_name, blk_dim in solver.gf_struct:
        Sigma_HartreeFock[blk_name] = np.diag(hf_per_color[c:c + blk_dim]).astype(float)
        c += blk_dim
    return Sigma_HartreeFock


def _legendre_discontinuity_weights(g_l):
    beta = g_l.mesh.beta
    l = np.arange(g_l.data.shape[0])
    weights = np.zeros_like(l, dtype=float)
    even_l = (l % 2) == 0
    weights[even_l] = -2.0 * np.sqrt(2.0 * l[even_l] + 1.0) / beta
    return weights


def _legendre_discontinuity(g_l):
    weights = _legendre_discontinuity_weights(g_l)
    return np.tensordot(weights, g_l.data, axes=(0, 0))


def _enforce_legendre_discontinuity(g_l, discontinuity):
    """Complex-valued equivalent of GfLegendre.enforce_discontinuity."""
    weights = _legendre_discontinuity_weights(g_l)
    norm = np.dot(weights, weights)
    correction = np.asarray(discontinuity, dtype=complex) - _legendre_discontinuity(g_l)
    g_l.data[:] += correction[None, :, :] * (weights / norm)[:, None, None]


# =============================================================================
# Self-energy post-processing
# =============================================================================

def postprocess_sigma(
    solver,
    symmetrize_func=None,
    **post_proc_params
):
    """
    Post-process the self-energy from the CT-SEG solver.

    Steps:
    1. Fourier transform G(tau) to G(iw)
    2. Compute fermionic Weiss field g(iw) from Delta(tau)
    3. Compute static self-energy Sigma_HF (analytically or via tail fitting)
    4. Compute dynamic self-energy Sigma(iw) (Dyson or improved estimator)
    5. Tail fitting on Sigma(iw)
    6. Extract dynamic part (Sigma - Sigma_HF)
    """
    from triqs.gfs import iOmega_n
    from triqs.gfs.tools import inverse

    degenerate_blk = post_proc_params['degenerate_blk']

    def symmetrize(gf):
        if symmetrize_func is not None and degenerate_blk:
            gf << symmetrize_func(gf, degenerate_blk)

    mesh = MeshImFreq(beta=solver.beta, statistic="Fermion", n_iw=solver.n_iw)
    Sigma_iw = BlockGf(mesh=mesh, gf_struct=solver.gf_struct)
    Sigma_iw.zero()
    G_iw = Sigma_iw.copy()
    G0_iw = Sigma_iw.copy()
    F_iw = None

    # 1. Fourier transform G(tau) to G(iw)
    Gf_known_moments = make_zero_tail(G_iw, n_moments=2)
    for i, bl in enumerate(G_iw.indices):
        Gf_known_moments[i][1] = np.eye(G_iw[bl].target_shape[0])
        G_iw[bl].set_from_fourier(solver.results.G_tau[bl], Gf_known_moments[i])
    G_iw << make_hermitian(G_iw)
    symmetrize(G_iw)
    G_iw_from_tau = G_iw.copy()

    # 2. Compute fermionic Weiss field g(iw)
    Delta_iw = BlockGf(mesh=mesh, gf_struct=solver.gf_struct)
    Delta_known_moments = make_zero_tail(Delta_iw, n_moments=1)
    for i, bl in enumerate(solver.Delta_tau.indices):
        Delta_iw[bl].set_from_fourier(solver.Delta_tau[bl], Delta_known_moments[i])
        G0_iw[bl] << inverse(iOmega_n - solver.h_loc0_mat[i] - Delta_iw[bl])
    G0_iw << make_hermitian(G0_iw)
    symmetrize(G0_iw)

    def _report_sigma_hf(hf_dict):
        for blk_name, hf_val in hf_dict.items():
            mpi.report(f"Sigma_HF {blk_name}:")
            mpi.report(f"    {hf_val}")
        mpi.report("")

    # 3. Compute the HF self-energy
    if post_proc_params['analytic_hf']:
        Sigma_HartreeFock = compute_sigma_hartreefock(solver)

        if symmetrize_func is not None and degenerate_blk:
            Sigma_HF_list = symmetrize_func(list(Sigma_HartreeFock.values()), degenerate_blk)
            Sigma_HartreeFock = dict(zip(Sigma_HartreeFock.keys(), Sigma_HF_list))

        _report_sigma_hf(Sigma_HartreeFock)
        solver.Sigma_moments = {
            blk_name: np.array([hf_val], dtype=complex)
            for blk_name, hf_val in Sigma_HartreeFock.items()
        }
    else:
        Sigma_HartreeFock = None
        solver.Sigma_moments = None

    # 4. Compute self-energy
    if solver.results.F_tau is None:
        mpi.report("F(tau) is not measured -> Compute the self-energy via the Dyson equation.\n")
        Sigma_iw_dyson = inverse(G0_iw) - inverse(G_iw)
        Sigma_iw = Sigma_iw_dyson.copy()
    else:
        mpi.report("F(tau) is measured -> Compute the self-energy via the improved estimator.\n")
        Sigma_iw_dyson = inverse(G0_iw) - inverse(G_iw)
        F_iw = G_iw.copy()
        F_iw << 0.0
        F_known_moments = make_zero_tail(F_iw, n_moments=1)
        for i, bl in enumerate(F_iw.indices):
            F_iw[bl].set_from_fourier(solver.results.F_tau[bl], F_known_moments[i])
        F_iw << make_hermitian(F_iw)
        symmetrize(F_iw)

        for block, fw in F_iw:
            for iw in fw.mesh:
                Sigma_iw[block][iw] = fw[iw] / G_iw[block][iw]

    Sigma_iw << make_hermitian(Sigma_iw)
    symmetrize(Sigma_iw)

    if post_proc_params['perform_tail_fit']:
        Sigma_iw = tail_fit(
            Sigma_iw,
            fit_min_n=post_proc_params['fit_min_n'],
            fit_max_n=post_proc_params['fit_max_n'],
            fit_min_w=post_proc_params['fit_min_w'],
            fit_max_w=post_proc_params['fit_max_w'],
            fit_max_moment=post_proc_params['fit_max_moment'],
            fit_known_moments=None,
        )

    Sigma_iw << make_hermitian(Sigma_iw)
    symmetrize(Sigma_iw)

    # 5. Update G(iw) with the fitted self-energy
    G_iw << inverse(inverse(G0_iw) - Sigma_iw)
    G_iw << make_hermitian(G_iw)
    symmetrize(G_iw)

    G_l = G_iw_l = Sigma_iw_dyson_l = None
    if getattr(solver.results, 'G_l', None) is not None:
        G_l = solver.results.G_l.copy()
        G_iw_l = G_iw.copy()
        G_iw_l << 0.0
        for bl, g_l in G_l:
            _enforce_legendre_discontinuity(g_l, np.eye(g_l.target_shape[0]))
            G_iw_l[bl].set_from_legendre(g_l)
        G_iw_l << make_hermitian(G_iw_l)
        symmetrize(G_iw_l)
        Sigma_iw_dyson_l = inverse(G0_iw) - inverse(G_iw_l)
        Sigma_iw_dyson_l << make_hermitian(Sigma_iw_dyson_l)
        symmetrize(Sigma_iw_dyson_l)

    F_l_raw = F_l = F_iw_l = Sigma_iw_l = None
    if getattr(solver.results, 'F_l', None) is not None:
        # Unlike G_l, F_l has no universal discontinuity.  On this standalone
        # branch preserve the measured coefficients instead of depending on
        # analytic tail moments supplied by another feature branch.
        F_l_raw = solver.results.F_l
        F_l = solver.results.F_l.copy()
        F_iw_l = G_iw.copy()
        F_iw_l << 0.0
        for bl, f_l in F_l:
            F_iw_l[bl].set_from_legendre(f_l)
        F_iw_l << make_hermitian(F_iw_l)
        symmetrize(F_iw_l)

        G_for_F_l = G_iw_l if G_iw_l is not None else G_iw_from_tau
        Sigma_iw_l = Sigma_iw.copy()
        Sigma_iw_l << 0.0
        for bl, f_iw_l in F_iw_l:
            Sigma_iw_l[bl] << f_iw_l * inverse(G_for_F_l[bl])
        Sigma_iw_l << make_hermitian(Sigma_iw_l)
        symmetrize(Sigma_iw_l)

    if Sigma_HartreeFock is None:
        Sigma_HartreeFock = {block: gf.fit_hermitian_tail()[0][0] for block, gf in Sigma_iw}
        mpi.report("Extracting the static self-energy via tail fitting:")
        _report_sigma_hf(Sigma_HartreeFock)

    # 6. Extract the dynamic part (Sigma - Sigma_HF)
    Sigma_dynamic = Sigma_iw.copy()
    for bl, g in Sigma_dynamic:
        g << g - Sigma_HartreeFock[bl]

    return {
        'G_iw': G_iw,
        'G_iw_from_tau': G_iw_from_tau,
        'G_tau': solver.results.G_tau,
        'G_l': G_l,
        'G_iw_l': G_iw_l,
        'F_iw': F_iw,
        'F_l_raw': F_l_raw,
        'F_l': F_l,
        'F_iw_l': F_iw_l,
        'Sigma_iw': Sigma_iw,
        'Sigma_iw_dyson': Sigma_iw_dyson,
        'Sigma_iw_dyson_l': Sigma_iw_dyson_l,
        'Sigma_iw_l': Sigma_iw_l,
        'Sigma_dynamic': Sigma_dynamic,
        'Sigma_HartreeFock': list(Sigma_HartreeFock.values()),
    }


# =============================================================================
# Polarizability and screened interaction
# =============================================================================

def postprocess_pi(
    solver,
    degenerate_blk=None,
    symmetrize_func=None,
    output_in_4idx=False,
    truncate_uchi=False,
):
    """Post-process the charge susceptibility to obtain the impurity polarizability and screened interaction."""
    mpi.report("Charge susceptibility is measured for a impurity with dynamic interactions")
    mpi.report('--> Post-processing the density-density susceptibility to obtain the impurity polarizability.\n')

    def symmetrize(gf):
        if symmetrize_func is not None and degenerate_blk is not None:
            gf << symmetrize_func(gf, degenerate_blk)

    block_name, index_in_block, color_to_orbital = build_color_tables(solver.gf_struct)
    n_color = len(block_name)
    assert n_color % 2 == 0, "n_color is an odd number."
    n_orb = n_color // 2

    D0_tau = _assemble_color_gf(solver.D0_tau, block_name, index_in_block)
    nn_iw_dlr = _assemble_color_gf(solver.results.nn_nu_dlr, block_name, index_in_block)
    iw_mesh_dlr = nn_iw_dlr.mesh

    # 1. Subtract the constant part of the charge susceptibility
    densities = np.array([
        solver.results.densities[block_name[c]][index_in_block[c]]
        for c in range(n_color)
    ], dtype=float)
    mpi.report(f"Average of time-dependent occupations: {densities}")

    mpi.report("Subtracting the constant component, and then symmetrizing the density-density susceptibility: \n"
               "  1. nn(t).imag = 0.0\n"
               "  2. nn(i, j) = nn(j, i)\n")

    w0_idx = next(k for k, iw in enumerate(iw_mesh_dlr) if iw.index == 0)
    beta = iw_mesh_dlr.beta
    for c1 in range(n_color):
        for c2 in range(c1 + 1):
            dd = beta * densities[c1] * densities[c2]
            nn_iw_dlr[c1, c2].data[w0_idx] -= dd
            nn_iw_dlr[c1, c2].data.imag = 0.0
            if c1 != c2:
                nn_iw_dlr[c2, c1].data[w0_idx] -= dd
                nn_iw_dlr[c1, c2].data[:] += nn_iw_dlr[c2, c1].data[:].real
                nn_iw_dlr[c1, c2].data[:] /= 2.0
                nn_iw_dlr[c2, c1] << nn_iw_dlr[c1, c2]

    nn_iw = make_gf_imfreq(make_gf_dlr(nn_iw_dlr), n_iw=solver.n_iw)

    # Convert to density-density basis
    nn_iw_dd = Gf(mesh=nn_iw.mesh, target_shape=[n_orb, n_orb])
    for c1, c2 in product(range(n_color), repeat=2):
        nn_iw_dd[color_to_orbital[c1], color_to_orbital[c2]].data[:] += nn_iw[c1, c2].data[:]

    symmetrize(nn_iw_dd)

    # Convert to product basis
    nn_iw_pb = Gf(mesh=nn_iw.mesh, target_shape=[n_orb*n_orb, n_orb*n_orb])
    for i, j in product(range(n_orb), repeat=2):
        if i < j:
            continue
        ii, jj = i*n_orb+i, j*n_orb+j
        nn_iw_pb[ii, jj] << nn_iw_dd[i, j].real
        if i != j:
            nn_iw_pb[jj, ii] << nn_iw_pb[ii, jj]

    # 2. Construct bosonic Weiss field in product basis
    D0_iw = nn_iw.copy()
    D0_iw << 0.0
    u_known_moments = make_zero_tail(D0_iw, n_moments=2)
    D0_iw.set_from_fourier(D0_tau, u_known_moments)

    D0_iijj = D0_iw[:n_orb, n_orb:2*n_orb]
    D0_ijij = D0_iijj - D0_iw[:n_orb, 0:n_orb]

    Vijkl = extract_u_tensor_from_h_int(h_int=solver.h_int, gf_struct=solver.gf_struct, return_4idx=True)

    U_iw_pb = Gf(mesh=nn_iw_pb.mesh, target_shape=nn_iw_pb.target_shape)
    for i, j in product(range(n_orb), repeat=2):
        if i < j:
            continue
        ii, jj = i*n_orb+i, j*n_orb+j
        U_iw_pb[ii, jj] << D0_iijj[i, j].real
        U_iw_pb[ii, jj].data[:] += Vijkl[i, j, i, j]
        if i == j:
            continue
        ij, ji = i*n_orb+j, j*n_orb+i
        U_iw_pb[jj, ii] << U_iw_pb[ii, jj]
        U_iw_pb[ij, ij] << D0_ijij[i, j].real
        U_iw_pb[ij, ij].data[:] += Vijkl[i, j, j, i]
        U_iw_pb[ji, ji] << U_iw_pb[ij, ij]
        U_iw_pb[ij, ji] << D0_ijij[i, j].real
        U_iw_pb[ij, ji].data[:] += Vijkl[i, i, j, j]
        U_iw_pb[ji, ij] << U_iw_pb[ij, ji]

    # 3. Dyson equation: Pi(w) = Chi(w) * [U(w)*Chi(w) - I]^-1
    Pi_iw_pb = Gf(mesh=nn_iw_pb.mesh, target_shape=nn_iw_pb.target_shape)
    ones = np.eye(n_orb*n_orb, dtype=complex)
    for iwn in nn_iw_pb.mesh:
        UX = U_iw_pb[iwn] @ nn_iw_pb[iwn]
        UX = check_spectrum(UX, truncation=truncate_uchi)
        denom = UX - ones
        cond = np.linalg.cond(denom)
        if cond > 20:
            mpi.report(f"WARNING: Large condition number for [U(w) * Chi(w) - I] = {cond} at n = {iwn.index}.")
        Pi_iw_pb[iwn] = nn_iw_pb[iwn] @ np.linalg.pinv(denom)
        Pi_iw_pb[iwn].imag = 0.0

    # 4. Screened interaction W(w) = U(w) - U(w) * Chi(w) * U(w)
    W_iw_pb = Gf(mesh=nn_iw_pb.mesh, target_shape=nn_iw_pb.target_shape)
    for iwn in nn_iw_pb.mesh:
        W_iw_pb[iwn] = U_iw_pb[iwn] - U_iw_pb[iwn] @ nn_iw_pb[iwn] @ U_iw_pb[iwn]

    # 5. Remove the static part
    for i, j in product(range(n_orb), repeat=2):
        if i < j:
            continue
        ii, jj = i*n_orb+i, j*n_orb+j
        W_iw_pb[ii, jj].data[:] -= Vijkl[i, j, i, j]
        if i == j:
            continue
        ij, ji = i*n_orb+j, j*n_orb+i
        W_iw_pb[jj, ii] << W_iw_pb[ii, jj]
        W_iw_pb[ji, ji].data[:] -= Vijkl[i, j, j, i]
        W_iw_pb[ij, ij] << W_iw_pb[ji, ji]
        W_iw_pb[ji, ij].data[:] -= Vijkl[i, i, j, j]
        W_iw_pb[ij, ji] << W_iw_pb[ji, ij]

    # Transform back
    if not output_in_4idx:
        Chi_iw = nn_iw_dd
        Pi_iw = Gf(mesh=nn_iw_pb.mesh, target_shape=(n_orb, n_orb))
        W_iw = Pi_iw.copy()
        for i, j in product(range(n_orb), repeat=2):
            Pi_iw[i, j] << Pi_iw_pb[i*n_orb+i, j*n_orb+j]
            W_iw[i, j] << W_iw_pb[i*n_orb+i, j*n_orb+j]
        symmetrize(Pi_iw)
        symmetrize(W_iw)
    else:
        Chi_iw = nn_iw_pb
        Pi_iw = Gf(mesh=nn_iw_pb.mesh, target_shape=(n_orb, n_orb, n_orb, n_orb))
        W_iw = Pi_iw.copy()
        for i, j, k, l in product(range(n_orb), repeat=4):
            Pi_iw[i, j, k, l] << Pi_iw_pb[i*n_orb+j, k*n_orb+l]
            W_iw[i, j, k, l] << W_iw_pb[i*n_orb+j, k*n_orb+l]

    return {
        'Pi_iw': Pi_iw,
        'W_iw': W_iw,
        'Chi_iw': Chi_iw,
    }


# =============================================================================
# Top-level post-processing dispatcher
# =============================================================================

def postprocess(
    solver,
    symmetrize_func=None,
    **post_proc_params,
):
    pp_results = postprocess_sigma(solver, symmetrize_func, **post_proc_params)

    has_chi = solver.D0_tau is not None and getattr(solver.results, 'nn_nu_dlr', None) is not None
    if has_chi:
        deg_blk = post_proc_params['degenerate_blk']
        deg_blk_2e = (
            [np.array(blks[:len(blks) // 2]) for blks in deg_blk]
            if deg_blk else None
        )
        pp_results.update(
            postprocess_pi(
                solver, deg_blk_2e, symmetrize_func,
                truncate_uchi=post_proc_params['truncate_uchi'],
            )
        )

    wmax = post_proc_params.get('dlr_w_max')
    eps = post_proc_params.get('dlr_eps')
    if wmax is not None and eps is not None:
        from triqs.solver_utils import make_gf_dlr_imfreq
        keys = ['Sigma_iw', 'G_iw'] + (['Pi_iw', 'W_iw', 'Chi_iw'] if has_chi else [])
        for k in keys:
            pp_results[k] = make_gf_dlr_imfreq(pp_results[k], w_max=wmax, eps=eps)

    return pp_results
