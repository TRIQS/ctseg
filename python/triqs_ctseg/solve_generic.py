"""
Generic solver interface for CT-SEG.

Provides a functional API to the triqs_ctseg solver with:
- Dynamic dispatch based on mesh type (MeshImFreq, MeshDLRImFreq)
- Dynamic or static Coulomb interactions
- Optional post-processing for self-energy and polarization extraction
"""

import numpy as np
import triqs.utility.mpi as mpi

from triqs.gfs import (
    MeshDLRImFreq, MeshImFreq, make_gf_from_fourier,
    make_gf_imtime, fit_hermitian_tail
)
from triqs.gfs.tools import make_zero_tail
from triqs.operators import Operator
from triqs.operators.util.extractors import block_matrix_from_op, op_from_block_matrix
from triqs.solver_utils import SolverResults

from triqs_ctseg import Solver
from .postprocessing import postprocess


def _canonicalize_h_loc0(h_loc0, gf_struct):
    """Normalize h_loc0 into (Operator, block-matrix array) regardless of input form.

    Accepts either a many_body_operator or an iterable of dense block matrices
    (one per block in gf_struct).
    """
    if isinstance(h_loc0, Operator):
        return h_loc0, block_matrix_from_op(h_loc0, gf_struct)
    bl = np.empty(len(h_loc0), dtype=object)
    bl[:] = [np.asarray(m) for m in h_loc0]
    return op_from_block_matrix(bl, gf_struct), bl


_POST_PROC_DEFAULTS = {
    'perform_tail_fit': False,
    'fit_max_moment': 3,
    'fit_min_w': None,
    'fit_max_w': None,
    'fit_min_n': None,
    'fit_max_n': None,
    'analytic_hf': False,
    'degenerate_blk': None,
    'truncate_uchi': False,
}


def _prepare_solver(
    Delta_iw,
    h_loc0_bl,
    h_int,
    solve_density_only=False,
    **solver_interface_params,
):
    """Prepare the CT-SEG solver instance and split parameters into solver/post-processing."""
    gf_struct = [(bl, gf.target_shape[0]) for (bl, gf) in Delta_iw]
    mesh = Delta_iw.mesh
    beta = mesh.beta

    solver_params = solver_interface_params.copy()

    n_iw = solver_params.pop('n_iw', 1025)
    n_tau = solver_params.pop('n_tau', 10001)
    n_tau_bosonic = solver_params.pop('n_tau_bosonic', n_tau)
    n_l = solver_params.pop('n_l', None)

    solver_params['measure_densities'] = True
    solver_params['measure_G_tau'] = not solve_density_only
    solver_params['measure_F_tau'] = not solve_density_only
    solver_params['measure_nn_tau'] = not solve_density_only

    if isinstance(mesh, MeshDLRImFreq):
        solver_params.setdefault('dlr_omega_max', mesh.w_max)
        solver_params.setdefault('dlr_eps', mesh.eps)
        max_dlr_idx = max(abs(iw.index) for iw in mesh)
        if max_dlr_idx > n_iw:
            mpi.report(f"WARNING: n_iw = {n_iw} is smaller than the maximum DLR frequency index"
                       f" ({max_dlr_idx}). Setting n_iw to {max_dlr_idx+1}.")
            n_iw = max_dlr_idx + 1

    solver_constructor_params = {
        'gf_struct': gf_struct,
        'beta': beta,
        'n_tau': n_tau,
        'n_tau_bosonic': n_tau_bosonic,
    }
    if solver_params.get('measure_G_l', False) or solver_params.get('measure_F_l', False):
        solver_constructor_params['n_l'] = 30 if n_l is None else n_l
    S = Solver(**solver_constructor_params)
    # Solver does not expose these as attributes; postprocess_sigma relies on them
    S.n_iw = n_iw
    S.beta = beta
    S.gf_struct = gf_struct
    S.h_int = h_int
    S.h_loc0_mat = h_loc0_bl

    post_proc_params = {k: solver_params.pop(k, default) for k, default in _POST_PROC_DEFAULTS.items()}
    post_proc_params['dlr_w_max'] = solver_params.pop('dlr_omega_max', None)
    post_proc_params['dlr_eps'] = solver_params.pop('dlr_eps', None)

    return S, n_tau, solver_params, post_proc_params


def _prepare_delta_tau(S, Delta_iw, D0_iw, n_tau):
    """Prepare Delta(tau) and D0(tau) on the solver's tau mesh."""
    mesh = Delta_iw.mesh

    if isinstance(mesh, MeshDLRImFreq):
        S.Delta_tau << make_gf_imtime(Delta_iw, n_tau)
        if D0_iw is not None:
            for name1, name2 in D0_iw.indices:
                n_tau_b = len(S.D0_tau[name1, name2].mesh)
                S.D0_tau[name1, name2] << make_gf_imtime(D0_iw[name1, name2], n_tau_b)
    elif isinstance(mesh, MeshImFreq):
        for block, _ in S.Delta_tau:
            S.Delta_tau[block] << make_gf_from_fourier(
                Delta_iw[block],
                S.Delta_tau[block].mesh,
                fit_hermitian_tail(Delta_iw[block], make_zero_tail(Delta_iw[block], 1))[0],
            )
        if D0_iw is not None:
            for name1, name2 in D0_iw.indices:
                S.D0_tau[name1, name2] << make_gf_from_fourier(
                    D0_iw[name1, name2],
                    S.D0_tau[name1, name2].mesh,
                    fit_hermitian_tail(D0_iw[name1, name2], make_zero_tail(D0_iw[name1, name2], n_moments=2))[0],
                )
    else:
        raise NotImplementedError(f"Unsupported mesh type: {type(mesh)}")


# =============================================================================
# Main solve functions
# =============================================================================

def _run_solver(Delta_iw, h_loc0, h_int, D0_iw, solve_density_only, **params):
    """Set up the solver, run MC, and return (S, post_proc_params)."""
    gf_struct = [(bl, gf.target_shape[0]) for (bl, gf) in Delta_iw]
    h_loc0_op, h_loc0_bl = _canonicalize_h_loc0(h_loc0, gf_struct)

    S, n_tau, solver_params, post_proc_params = _prepare_solver(
        Delta_iw, h_loc0_bl, h_int,
        solve_density_only=solve_density_only, **params,
    )
    _prepare_delta_tau(S, Delta_iw, D0_iw, n_tau)

    mpi.report("Solving the impurity problem with CT-SEG"
               + (" for density" if solve_density_only else ""))
    S.solve(h_loc0=h_loc0_op, h_int=h_int, **solver_params)
    return S, post_proc_params


def solve_generic(
    Delta_iw,
    h_loc0,
    h_int,
    D0_iw=None,
    symmetrize_func=None,
    **solver_interface_params,
):
    """Solve the impurity problem using CT-SEG.

    Parameters
    ----------
    Delta_iw : BlockGf
        Hybridization function on MeshImFreq or MeshDLRImFreq.
    h_loc0 : Operator | list[np.ndarray]
        Local non-interacting Hamiltonian. Accepted as either a many_body_operator
        or an iterable of dense block matrices (one per block in Delta_iw).
    h_int : Operator
        Interaction Hamiltonian.
    D0_iw : Block2Gf, optional
        Bosonic Weiss field for dynamic interactions.
    symmetrize_func : callable, optional
        Symmetrization function for post-processing.
    **solver_interface_params
        Solver and post-processing parameters (n_iw, n_tau, n_cycles, etc.).

    Returns
    -------
    SolverResults
    """
    S, post_proc_params = _run_solver(
        Delta_iw, h_loc0, h_int, D0_iw,
        solve_density_only=False, **solver_interface_params,
    )

    pp_results = postprocess(S, symmetrize_func=symmetrize_func, **post_proc_params)
    known_result_fields = set(SolverResults.__dataclass_fields__)
    result_kwargs = {
        k: v for k, v in pp_results.items()
        if v is not None and k in known_result_fields
    }
    extra_kwargs = {
        k: v for k, v in pp_results.items()
        if v is not None and k not in known_result_fields
    }
    result_kwargs['Solver'] = S
    results = SolverResults(**result_kwargs)
    for k, v in extra_kwargs.items():
        setattr(results, k, v)
    return results


def solve_density(
    Delta_iw,
    h_loc0,
    h_int,
    D0_iw=None,
    **solver_interface_params,
):
    """Solve the impurity problem with CT-SEG, measuring only densities.

    Parameters
    ----------
    Delta_iw : BlockGf
        Hybridization function on MeshImFreq or MeshDLRImFreq.
    h_loc0 : Operator | list[np.ndarray]
        Local non-interacting Hamiltonian. Accepted as either a many_body_operator
        or an iterable of dense block matrices (one per block in Delta_iw).
    h_int : Operator
        Interaction Hamiltonian.
    D0_iw : Block2Gf, optional
        Bosonic Weiss field for dynamic interactions.
    **solver_interface_params
        Solver parameters (n_iw, n_tau, n_cycles, etc.).

    Returns
    -------
    SolverResults
    """
    S, _ = _run_solver(
        Delta_iw, h_loc0, h_int, D0_iw,
        solve_density_only=True, **solver_interface_params,
    )
    return SolverResults(Solver=S)
