# Copyright (c) 2026, The Simons Foundation
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

"""MeshDLRImFreq + dynamic (retarded) interaction end-to-end test.

This combination is the only code path in solve_generic that is not
already covered by solve_generic_dlr.py or solve_generic_dynamical_U.py:
the `if D0_iw is not None:` branch inside the MeshDLRImFreq block of
_prepare_delta_tau, which funnels the bosonic Weiss field through
make_gf_imtime instead of Fourier + tail fit.
"""

from triqs.gf import (
    MeshDLRImFreq, Gf, BlockGf, Block2Gf, inverse, iOmega_n,
)
from triqs.gf.descriptors import Function
from triqs.operators import n
import triqs.utility.mpi as mpi
import h5
from triqs.utility.h5diff import h5diff

from triqs_ctseg import solve_generic

beta = 10.
U = 4.
mu = U / 2.
eps = 0.2
wp = 1.
g_coupling = 1.
w_max = 10.
dlr_eps = 1e-10
n_tau = 2051
n_tau_bosonic = 2001

# Fermionic hybridization on MeshDLRImFreq
mesh_f = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=dlr_eps, symmetrize=True)
gf_struct = [('down', 1), ('up', 1)]
Delta_iw = BlockGf(
    name_block_generator=[(bl, Gf(mesh=mesh_f, target_shape=[dim, dim]))
                          for bl, dim in gf_struct],
    make_copies=False,
)
Delta_iw << inverse(iOmega_n - eps)

# Bosonic Weiss field D0_iw as a Block2Gf on MeshDLRImFreq
mesh_b = MeshDLRImFreq(beta=beta, statistic='Boson', w_max=w_max, eps=dlr_eps, symmetrize=True)
block_names = [bl for bl, _ in gf_struct]
D0_blocks = [[Gf(mesh=mesh_b, target_shape=[1, 1]) for _ in block_names]
             for _ in block_names]
D0_iw = Block2Gf(block_names, block_names, D0_blocks, make_copies=False)
for n1, n2 in D0_iw.indices:
    D0_iw[n1, n2] << Function(lambda w: g_coupling * wp**2 / (w**2 - wp**2))

h_loc0 = -mu * (n('up', 0) + n('down', 0))
h_int = U * n('up', 0) * n('down', 0)

results = solve_generic(
    Delta_iw, h_loc0, h_int, D0_iw=D0_iw,
    n_tau=n_tau,
    n_tau_bosonic=n_tau_bosonic,
    length_cycle=50,
    n_warmup_cycles=1000,
    n_cycles=10000,
    analytic_hf=True,
)

if mpi.is_master_node():
    # DLR output conversion should have been applied
    sigma_mesh_type = type(results.Sigma_iw.mesh).__name__
    assert 'DLR' in sigma_mesh_type, f"expected DLR mesh, got {sigma_mesh_type}"

    # Analytic HF + dynamic D0 exercises extract_screen_matrix_from_D0_tau
    assert results.Sigma_HartreeFock is not None

    with h5.HDFArchive("solve_generic_dlr_dynamical_U.out.h5", 'w') as A:
        A['G_tau'] = results.G_tau
        A['G_iw'] = results.G_iw
        A['Sigma_iw'] = results.Sigma_iw
        A['Sigma_dynamic'] = results.Sigma_dynamic
        A['Sigma_HartreeFock'] = results.Sigma_HartreeFock
        A['densities'] = results.Solver.results.densities

    h5diff("solve_generic_dlr_dynamical_U.out.h5",
           "solve_generic_dlr_dynamical_U.ref.h5", precision=1e-9)
