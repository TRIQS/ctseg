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


def _block2gf_is_nonzero(gf, tol=1e-13):
    if gf is None:
        return False
    return any(np.max(np.abs(gf[key].data)) > tol for key in gf.indices)


def _gf_is_nonzero(gf, tol=1e-13):
    if gf is None:
        return False
    return np.max(np.abs(gf.data)) > tol


def _density_observables_from_density_matrix(solver):
    hist = getattr(solver, 'density_matrix', None)
    if hist is None:
        hist = getattr(solver.results, 'state_hist', None)
    if hist is None:
        return None, None

    prob = np.asarray(hist, dtype=float)
    block_name, _, _ = build_color_tables(solver.gf_struct)
    n_color = len(block_name)
    if prob.size != 2 ** n_color:
        raise ValueError(
            f"density matrix has length {prob.size}, expected {2 ** n_color} for {n_color} colors"
        )

    states = np.arange(prob.size, dtype=np.int64)
    occ = ((states[:, None] >> np.arange(n_color)) & 1).astype(float)
    n_avg = prob @ occ
    nn_avg = np.einsum('s,sc,sd->cd', prob, occ, occ)
    return n_avg, nn_avg


def _densities_from_results(solver):
    densities = getattr(solver.results, 'densities', None)
    if densities is None:
        return None
    block_name, index_in_block, _ = build_color_tables(solver.gf_struct)
    return np.array([
        densities[block_name[c]][index_in_block[c]]
        for c in range(len(block_name))
    ], dtype=float)


def _nn_static_from_results(solver):
    nn_static = getattr(solver.results, 'nn_static', None)
    if nn_static is None:
        return None

    block_name, index_in_block, _ = build_color_tables(solver.gf_struct)
    n_color = len(block_name)
    out = np.zeros((n_color, n_color), dtype=float)
    for c1, c2 in product(range(n_color), repeat=2):
        out[c1, c2] = nn_static[block_name[c1], block_name[c2]][index_in_block[c1], index_in_block[c2]]
    return out


def _color_vector_to_block_diag(vec, gf_struct):
    out = {}
    offset = 0
    for blk_name, blk_dim in gf_struct:
        out[blk_name] = np.diag(vec[offset:offset + blk_dim]).astype(complex)
        offset += blk_dim
    return out


def _tail_dict_from_block_moments(block_moments):
    return {
        blk_name: np.asarray(moments, dtype=complex)
        for blk_name, moments in block_moments.items()
    }


def _d0_tau0_matrix(solver):
    block_name, index_in_block, _ = build_color_tables(solver.gf_struct)
    D0_tau = _assemble_color_gf(solver.D0_tau, block_name, index_in_block)
    return 0.5 * (D0_tau.data[0].real + D0_tau.data[-1].real)


def _spin_color_indices(gf_struct):
    block_name, _, _ = build_color_tables(gf_struct)
    up = [
        c for c, name in enumerate(block_name)
        if name.lower().startswith('up')
    ]
    down = [
        c for c, name in enumerate(block_name)
        if name.lower().startswith('down') or name.lower().startswith('dn')
    ]
    if len(up) != 1 or len(down) != 1:
        return None
    return up[0], down[0]


def _h_loc0_color_matrix(solver):
    block_name, _, _ = build_color_tables(solver.gf_struct)
    h_loc0_color = np.zeros((len(block_name), len(block_name)), dtype=complex)
    offset = 0
    for h_block in solver.h_loc0_mat:
        dim = h_block.shape[0]
        h_loc0_color[offset:offset + dim, offset:offset + dim] = h_block
        offset += dim
    return h_loc0_color


def _is_unpolarized_single_orbital(solver, up, down, tol=1e-12):
    block_name, _, _ = build_color_tables(solver.gf_struct)
    if len(block_name) != 2:
        return False
    h_loc0 = _h_loc0_color_matrix(solver)
    return abs(h_loc0[up, up] - h_loc0[down, down]) < tol


def _trapezoid_uniform(values, beta):
    values = np.asarray(values)
    if values.size < 2:
        raise ValueError("Need at least two tau points for trapezoidal integration.")
    dtau = beta / (values.size - 1)
    weights = np.ones(values.size)
    weights[0] = weights[-1] = 0.5
    return dtau * np.sum(weights * values)


def _periodic_tau_convolution_zero(left, middle, right, beta):
    """Approximate int dt dt' left(t) middle(t-t') right(t') / beta.

    The bosonic tau meshes include both 0 and beta. For the circular
    convolution use the half-open grid [0, beta) to avoid double-counting the
    endpoint.
    """
    left = np.asarray(left[:-1], dtype=float)
    middle = np.asarray(middle[:-1], dtype=float)
    right = np.asarray(right[:-1], dtype=float)
    if not (left.size == middle.size == right.size):
        raise ValueError("Convolution inputs must live on the same tau mesh.")
    n_tau = left.size
    dtau = beta / n_tau
    total = 0.0
    for i in range(n_tau):
        total += left[i] * np.dot(middle[(i - np.arange(n_tau)) % n_tau], right)
    return (dtau * dtau / beta) * total


def _phase2c_jperp_components_from_tau(jperp_tau, chi_xx_tau, beta, U):
    jperp_tau = np.asarray(jperp_tau, dtype=float)
    chi_xx_tau = np.asarray(chi_xx_tau, dtype=float)
    if jperp_tau.shape != chi_xx_tau.shape:
        raise ValueError("Jperp_tau and chi_xx_tau must have the same tau mesh.")

    j_tau0 = 0.5 * (jperp_tau[0] + jperp_tau[-1])
    j_chi = _trapezoid_uniform(jperp_tau * chi_xx_tau, beta)
    j_chi_j = _periodic_tau_convolution_zero(jperp_tau, chi_xx_tau, jperp_tau, beta)

    pure = 0.5 * (j_tau0 + j_chi_j)
    mixed = -2.0 * U * j_chi
    return {
        'pure': pure,
        'mixed': mixed,
        'total': pure + mixed,
        'int_J_chi_xx': j_chi,
        'J_chi_xx_J': j_chi_j,
        'J_tau0': j_tau0,
    }


def _phase2c_jperp_oriented_components_from_tau(
    jperp_tau,
    chi_minus_plus_tau,
    chi_plus_minus_tau,
    beta,
    U,
):
    jperp_tau = np.asarray(jperp_tau, dtype=float)
    chi_minus_plus_tau = np.asarray(chi_minus_plus_tau, dtype=float)
    chi_plus_minus_tau = np.asarray(chi_plus_minus_tau, dtype=float)
    if jperp_tau.shape != chi_minus_plus_tau.shape or jperp_tau.shape != chi_plus_minus_tau.shape:
        raise ValueError("Jperp_tau and oriented transverse-spin correlators must have the same tau mesh.")

    j_tau0 = 0.5 * (jperp_tau[0] + jperp_tau[-1])

    def one_orientation(chi_tau):
        j_chi = _trapezoid_uniform(jperp_tau * chi_tau, beta)
        j_chi_j = _periodic_tau_convolution_zero(jperp_tau, chi_tau, jperp_tau, beta)
        pure = 0.25 * (2.0 * j_tau0 + j_chi_j)
        mixed = -U * j_chi
        return {
            'pure': pure,
            'mixed': mixed,
            'total': pure + mixed,
            'int_J_chi': j_chi,
            'J_chi_J': j_chi_j,
            'J_tau0': j_tau0,
        }

    up = one_orientation(chi_minus_plus_tau)
    down = one_orientation(chi_plus_minus_tau)
    return {
        'mode': 'oriented',
        'up': up,
        'down': down,
        'Sminus_Splus': up,
        'Splus_Sminus': down,
        'total_up': up['total'],
        'total_down': down['total'],
        'J_tau0': j_tau0,
    }


def _chi_xx_tau_from_solver(solver, up, down):
    sperp_tau = getattr(solver.results, 'Sperp_tau', None)
    if sperp_tau is not None:
        return np.asarray(sperp_tau.data[:, 0, 0].real)

    # A spin-unpolarized one-body Hamiltonian does not prove spin-rotational
    # invariance: anisotropic interactions can make chi_zz != chi_xx.  Without
    # a measured transverse correlator there is no safe reconstruction here.
    return None


def _oriented_sperp_tau_from_solver(solver):
    sm_sp_tau = getattr(solver.results, 'Sminus_Splus_tau', None)
    sp_sm_tau = getattr(solver.results, 'Splus_Sminus_tau', None)
    if sm_sp_tau is None or sp_sm_tau is None:
        return None
    return (
        np.asarray(sm_sp_tau.data[:, 0, 0].real),
        np.asarray(sp_sm_tau.data[:, 0, 0].real),
    )


def _assemble_jperp_phase2c_sigma1(solver, U):
    spin_colors = _spin_color_indices(solver.gf_struct)
    if spin_colors is None:
        mpi.report("WARNING: Jperp Phase 2c moments require a single up/down orbital; "
                   "skipping analytic transverse Sigma_1/F_2.")
        return None

    up, down = spin_colors
    is_unpolarized = _is_unpolarized_single_orbital(solver, up, down)

    jperp_tau = np.asarray(solver.Jperp_tau.data[:, 0, 0].real)
    U_spin = 0.5 * (U[up, down] + U[down, up])

    oriented_tau = _oriented_sperp_tau_from_solver(solver)
    if oriented_tau is not None:
        chi_minus_plus_tau, chi_plus_minus_tau = oriented_tau
        if jperp_tau.shape == chi_minus_plus_tau.shape and jperp_tau.shape == chi_plus_minus_tau.shape:
            components = _phase2c_jperp_oriented_components_from_tau(
                jperp_tau,
                chi_minus_plus_tau,
                chi_plus_minus_tau,
                solver.Jperp_tau.mesh.beta,
                U_spin,
            )

            correction = np.zeros(U.shape[0], dtype=float)
            correction[up] = components['up']['total']
            correction[down] = components['down']['total']
            return correction, components

        mpi.report("WARNING: Jperp_tau and asymmetric transverse-spin tau meshes differ; "
                   "skipping asymmetric transverse Sigma_1/F_2.")
        if not is_unpolarized:
            return None

    if not is_unpolarized:
        mpi.report("WARNING: spin-polarized Jperp moments require Sminus_Splus_tau and Splus_Sminus_tau; "
                   "skipping analytic transverse Sigma_1/F_2.")
        return None

    chi_xx_tau = _chi_xx_tau_from_solver(solver, up, down)
    if chi_xx_tau is None:
        mpi.report("WARNING: Jperp Phase 2c moments require nn_tau or Sperp_tau; "
                   "skipping analytic transverse Sigma_1/F_2.")
        return None

    if jperp_tau.shape != chi_xx_tau.shape:
        mpi.report("WARNING: Jperp_tau and chi_xx_tau tau meshes differ; "
                   "skipping analytic transverse Sigma_1/F_2.")
        return None

    components = _phase2c_jperp_components_from_tau(
        jperp_tau, chi_xx_tau, solver.Jperp_tau.mesh.beta, U_spin
    )

    correction = np.zeros(U.shape[0], dtype=float)
    correction[up] = components['total']
    correction[down] = components['total']
    return correction, components


def _assemble_density_tail_moments(solver, use_tail_moments=True):
    """Assemble diagonal color-space Sigma/G/F moments from static and D0 data."""
    if not use_tail_moments:
        return None

    n_avg, nn_avg = _density_observables_from_density_matrix(solver)
    if n_avg is None:
        n_avg = _densities_from_results(solver)
    if n_avg is None:
        mpi.report("WARNING: Cannot assemble tail moments because densities are not measured.")
        return None
    if nn_avg is None:
        nn_avg = _nn_static_from_results(solver)

    U = extract_u_tensor_from_h_int(h_int=solver.h_int, gf_struct=solver.gf_struct).real
    sigma0 = U @ n_avg
    sigma1 = None
    if nn_avg is not None:
        cov_nn = nn_avg - np.outer(n_avg, n_avg)
        sigma1 = np.einsum('cd,ce,de->c', U, U, cov_nn)

    has_D0 = _block2gf_is_nonzero(solver.D0_tau)
    if has_D0:
        D0_w0 = extract_screen_matrix_from_D0_tau(blk2_D0_tau=solver.D0_tau, gf_struct=solver.gf_struct).real
        phi_avg = D0_w0 @ n_avg
        sigma0 = sigma0 + phi_avg

        dyn_phi_n = getattr(solver.results, 'dyn_phi_n', None)
        dyn_phi_phi = getattr(solver.results, 'dyn_phi_phi', None)
        if dyn_phi_n is None or dyn_phi_phi is None:
            mpi.report("WARNING: D0_tau is non-zero but dyn_phi_n/dyn_phi_phi were not measured; "
                       "using Sigma_0/G_2 moments only and skipping analytic Sigma_1/F_2.")
            sigma1 = None
        elif nn_avg is None:
            mpi.report("WARNING: D0_tau is non-zero but density covariance is unavailable; "
                       "using Sigma_0/G_2 moments only and skipping analytic Sigma_1/F_2.")
            sigma1 = None
        else:
            dyn_phi_n = np.asarray(dyn_phi_n, dtype=float)
            dyn_phi_phi = np.asarray(dyn_phi_phi, dtype=float)
            cov_phi_n = dyn_phi_n - np.outer(phi_avg, n_avg)
            var_phi = np.diag(dyn_phi_phi) - phi_avg ** 2 + np.diag(_d0_tau0_matrix(solver))
            sigma1 = sigma1 + 2.0 * np.einsum('cd,cd->c', U, cov_phi_n) + var_phi

    has_Jperp = _gf_is_nonzero(getattr(solver, 'Jperp_tau', None))
    if has_Jperp:
        jperp_phase2c = _assemble_jperp_phase2c_sigma1(solver, U)
        if jperp_phase2c is not None and sigma1 is not None:
            jperp_sigma1, jperp_components = jperp_phase2c
            sigma1 = sigma1 + jperp_sigma1
            solver.Jperp_moment_components = jperp_components
        elif jperp_phase2c is not None:
            mpi.report("WARNING: static density covariance is unavailable; skipping analytic Jperp Sigma_1/F_2.")

    block_name, _, _ = build_color_tables(solver.gf_struct)
    h_loc0_color = _h_loc0_color_matrix(solver)

    sigma0_mat = np.diag(sigma0).astype(complex)
    g2_color = h_loc0_color + sigma0_mat

    Sigma_moments = {}
    G_moments = {}
    F_moments = {}
    F_tail_moments = {}
    offset = 0
    for blk_name, blk_dim in solver.gf_struct:
        sigma0_block = sigma0_mat[offset:offset + blk_dim, offset:offset + blk_dim]
        g2_block = g2_color[offset:offset + blk_dim, offset:offset + blk_dim]

        sigma_tail = np.zeros((1 if sigma1 is None else 2, blk_dim, blk_dim), dtype=complex)
        sigma_tail[0] = sigma0_block
        if sigma1 is not None:
            sigma_tail[1] = np.diag(sigma1[offset:offset + blk_dim]).astype(complex)
        Sigma_moments[blk_name] = sigma_tail

        g_tail = np.zeros((3, blk_dim, blk_dim), dtype=complex)
        g_tail[1] = np.eye(blk_dim)
        g_tail[2] = g2_block
        G_moments[blk_name] = g_tail

        F_moments[blk_name] = sigma0_block
        f_tail = np.zeros((2 if sigma1 is None else 3, blk_dim, blk_dim), dtype=complex)
        f_tail[1] = sigma0_block
        if sigma1 is not None:
            sigma1_block = np.diag(sigma1[offset:offset + blk_dim]).astype(complex)
            f_tail[2] = sigma0_block @ g2_block + sigma1_block
        F_tail_moments[blk_name] = f_tail
        offset += blk_dim

    return {
        'Sigma_HartreeFock': _color_vector_to_block_diag(sigma0, solver.gf_struct),
        'Sigma_moments': _tail_dict_from_block_moments(Sigma_moments),
        'G_moments': _tail_dict_from_block_moments(G_moments),
        'F_moments': F_moments,
        'F_tail_moments': _tail_dict_from_block_moments(F_tail_moments),
    }


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

    use_tail_moments = post_proc_params.get('use_tail_moments', True)
    tail_moments = _assemble_density_tail_moments(
        solver,
        use_tail_moments=use_tail_moments,
    )
    if tail_moments is not None:
        Sigma_HartreeFock = tail_moments['Sigma_HartreeFock']
        solver.Sigma_moments = tail_moments['Sigma_moments']
        solver.G_moments = tail_moments['G_moments']
        solver.F_moments = tail_moments['F_moments']
        solver.F_tail_moments = tail_moments['F_tail_moments']
    else:
        Sigma_HartreeFock = None
        solver.Sigma_moments = None
        solver.G_moments = None
        solver.F_moments = None
        solver.F_tail_moments = None

    mesh = MeshImFreq(beta=solver.beta, statistic="Fermion", n_iw=solver.n_iw)
    Sigma_iw = BlockGf(mesh=mesh, gf_struct=solver.gf_struct)
    Sigma_iw.zero()
    G_iw = Sigma_iw.copy()
    G0_iw = Sigma_iw.copy()

    # 1. Fourier transform G(tau) to G(iw)
    if not use_tail_moments or solver.G_moments is None:
        Gf_known_moments = make_zero_tail(G_iw, n_moments=2)
        for i, bl in enumerate(G_iw.indices):
            Gf_known_moments[i][1] = np.eye(G_iw[bl].target_shape[0])
            G_iw[bl].set_from_fourier(solver.results.G_tau[bl], Gf_known_moments[i])
    else:
        for bl in G_iw.indices:
            G_iw[bl].set_from_fourier(solver.results.G_tau[bl], solver.G_moments[bl])
    G_iw << make_hermitian(G_iw)
    symmetrize(G_iw)

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
    if Sigma_HartreeFock is None and post_proc_params['analytic_hf']:
        Sigma_HartreeFock = compute_sigma_hartreefock(solver)
        solver.Sigma_moments = {
            blk_name: np.array([hf_val], dtype=complex)
            for blk_name, hf_val in Sigma_HartreeFock.items()
        }

    if Sigma_HartreeFock is not None:
        if symmetrize_func is not None and degenerate_blk:
            Sigma_HF_list = symmetrize_func(list(Sigma_HartreeFock.values()), degenerate_blk)
            Sigma_HartreeFock = dict(zip(Sigma_HartreeFock.keys(), Sigma_HF_list))
        _report_sigma_hf(Sigma_HartreeFock)

    # 4. Compute self-energy
    if solver.results.F_tau is None:
        mpi.report("F(tau) is not measured -> Compute the self-energy via the Dyson equation.\n")
        Sigma_iw = inverse(G0_iw) - inverse(G_iw)
    else:
        mpi.report("F(tau) is measured -> Compute the self-energy via the improved estimator.\n")
        F_iw = G_iw.copy()
        F_iw << 0.0
        if not use_tail_moments:
            F_known_moments = make_zero_tail(F_iw, n_moments=1)
            for i, bl in enumerate(F_iw.indices):
                F_iw[bl].set_from_fourier(solver.results.F_tau[bl], F_known_moments[i])
        else:
            F_tau_for_fourier = solver.results.F_tau.copy()
            for bl, f_tau in F_tau_for_fourier:
                if solver.F_tail_moments is not None:
                    F_known_moments = solver.F_tail_moments[bl]
                else:
                    F_known_moments = np.zeros((2,) + f_tau.target_shape, dtype=complex)
                    F_known_moments[1] = -(f_tau.data[0] + f_tau.data[-1])
                f1 = F_known_moments[1]
                f_tau.data[0, :, :] = 0.5 * (f_tau.data[0, :, :] - f1 - f_tau.data[-1, :, :])
                f_tau.data[-1, :, :] = -f1 - f_tau.data[0, :, :]
                F_iw[bl].set_from_fourier(f_tau, F_known_moments)
        F_iw << make_hermitian(F_iw)
        symmetrize(F_iw)

        for block, fw in F_iw:
            if not use_tail_moments:
                for iw in fw.mesh:
                    Sigma_iw[block][iw] = fw[iw] / G_iw[block][iw]
            else:
                Sigma_iw[block] << fw * inverse(G_iw[block])

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
            fit_known_moments=solver.Sigma_moments if use_tail_moments else None,
        )

    Sigma_iw << make_hermitian(Sigma_iw)
    symmetrize(Sigma_iw)

    # 5. Update G(iw) with the fitted self-energy
    G_iw << inverse(inverse(G0_iw) - Sigma_iw)
    G_iw << make_hermitian(G_iw)
    symmetrize(G_iw)

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
        'G_tau': solver.results.G_tau,
        'Sigma_iw': Sigma_iw,
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
