# Copyright (c) 2026, The Simons Foundation
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

"""End-to-end test of solve_generic on MeshDLRImFreq with static interaction.

Exercises:
  - _prepare_delta_tau MeshDLRImFreq branch (make_gf_imtime, no tail fit)
  - compute_sigma_hartreefock via analytic_hf=True (uses
    extract_u_tensor_from_h_int + densities)
  - DLR output conversion in postprocess (dlr_w_max + dlr_eps pulled
    from the Delta_iw mesh -> make_gf_dlr_imfreq applied to
    Sigma_iw / G_iw)
"""

from triqs.gfs import MeshDLRImFreq, Gf, BlockGf, iOmega_n, inverse
from triqs.operators import n
import triqs.utility.mpi as mpi
import h5
from triqs.utility.h5diff import h5diff

from triqs_ctseg import solve_generic

# Model
beta = 10.
U = 4.
mu = U / 2.
eps = 0.2
w_max = 10.
dlr_eps = 1e-10
n_tau = 2051

# Hybridization on MeshDLRImFreq
dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=dlr_eps, symmetrize=True)
gf_struct = [('down', 1), ('up', 1)]
Delta_iw = BlockGf(
    name_block_generator=[(bl, Gf(mesh=dlr_iw_mesh, target_shape=[dim, dim]))
                          for bl, dim in gf_struct],
    make_copies=False,
)
Delta_iw << inverse(iOmega_n - eps)

h_loc0 = -mu * (n('up', 0) + n('down', 0))
h_int = U * n('up', 0) * n('down', 0)

results = solve_generic(
    Delta_iw, h_loc0, h_int,
    n_tau=n_tau,
    length_cycle=50,
    n_warmup_cycles=1000,
    n_cycles=10000,
    analytic_hf=True,
)

if mpi.is_master_node():
    # DLR conversion should have been applied: Sigma_iw / G_iw now live on
    # a DLR-imfreq mesh
    sigma_mesh_type = type(results.Sigma_iw.mesh).__name__
    g_mesh_type = type(results.G_iw.mesh).__name__
    assert 'DLR' in sigma_mesh_type, f"expected DLR mesh, got {sigma_mesh_type}"
    assert 'DLR' in g_mesh_type, f"expected DLR mesh, got {g_mesh_type}"

    # analytic_hf=True path populates solver.Sigma_moments with the HF values
    assert results.Solver.Sigma_moments is not None
    assert results.Sigma_HartreeFock is not None

    with h5.HDFArchive("solve_generic_dlr.out.h5", 'w') as A:
        A['G_tau'] = results.G_tau
        A['G_iw'] = results.G_iw
        A['Sigma_iw'] = results.Sigma_iw
        A['Sigma_dynamic'] = results.Sigma_dynamic
        A['Sigma_HartreeFock'] = results.Sigma_HartreeFock
        A['densities'] = results.Solver.results.densities

    h5diff("solve_generic_dlr.out.h5", "solve_generic_dlr.ref.h5", precision=1e-9)
