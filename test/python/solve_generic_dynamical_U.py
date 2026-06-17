# Copyright (c) 2026, The Simons Foundation
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

"""End-to-end test of solve_generic with a dynamic (retarded) interaction.

Mirrors the physics of dynamical_U.py.  Exercises:
  - _prepare_delta_tau D0_iw handling (Fourier + tail fit for MeshImFreq)
  - compute_sigma_hartreefock picking up the dynamic screening term via
    extract_screen_matrix_from_D0_tau (triggered by analytic_hf=True
    AND a non-None solver.D0_tau)

postprocess_pi does NOT run for ctseg because solver.results.nn_nu_dlr is
not measured; this test asserts that Pi/W/Chi are absent from the
returned SolverResults.
"""

from triqs.gfs import (
    MeshImFreq, Gf, BlockGf, Block2Gf, inverse, iOmega_n,
)
from triqs.gfs.descriptors import Function
from triqs.operators import n
import triqs.utility.mpi as mpi
import h5
from triqs.utility.h5diff import h5diff

from triqs_ctseg import solve_generic

# Model -- matches dynamical_U.py
beta = 10.
U = 4.
mu = U / 2.
eps = 0.2
wp = 1.
g = 1.
n_iw = 1025
n_tau = 2051
n_iw_bosonic = 1000
n_tau_bosonic = 2001

# Fermionic hybridization on MeshImFreq
iw_mesh_f = MeshImFreq(beta=beta, S='Fermion', n_iw=n_iw)
gf_struct = [('down', 1), ('up', 1)]
Delta_iw = BlockGf(
    name_block_generator=[(bl, Gf(mesh=iw_mesh_f, target_shape=[dim, dim]))
                          for bl, dim in gf_struct],
    make_copies=False,
)
Delta_iw << inverse(iOmega_n - eps)

# Bosonic Weiss field D0_iw as a Block2Gf on all (spin, spin') pairs
iw_mesh_b = MeshImFreq(beta=beta, S='Boson', n_iw=n_iw_bosonic)
block_names = [bl for bl, _ in gf_struct]
D0_blocks = [[Gf(mesh=iw_mesh_b, target_shape=[1, 1]) for _ in block_names]
             for _ in block_names]
D0_iw = Block2Gf(block_names, block_names, D0_blocks, make_copies=False)
for n1, n2 in D0_iw.indices:
    D0_iw[n1, n2] << Function(lambda w: g * wp**2 / (w**2 - wp**2))

h_loc0 = -mu * (n('up', 0) + n('down', 0))
h_int = U * n('up', 0) * n('down', 0)

results = solve_generic(
    Delta_iw, h_loc0, h_int, D0_iw=D0_iw,
    n_iw=n_iw,
    n_tau=n_tau,
    n_tau_bosonic=n_tau_bosonic,
    length_cycle=50,
    n_warmup_cycles=1000,
    n_cycles=10000,
    analytic_hf=True,
)

if mpi.is_master_node():
    # postprocess_pi should NOT have fired (nn_nu_dlr not measured in ctseg)
    assert 'Pi_iw' not in results.keys()
    assert 'W_iw' not in results.keys()
    assert 'Chi_iw' not in results.keys()

    # analytic_hf=True with dynamic interactions picks up the D0 screening
    # term in compute_sigma_hartreefock
    assert results.Sigma_HartreeFock is not None

    with h5.HDFArchive("solve_generic_dynamical_U.out.h5", 'w') as A:
        A['G_tau'] = results.G_tau
        A['G_iw'] = results.G_iw
        A['Sigma_iw'] = results.Sigma_iw
        A['Sigma_dynamic'] = results.Sigma_dynamic
        A['Sigma_HartreeFock'] = results.Sigma_HartreeFock
        A['densities'] = results.Solver.results.densities

    h5diff("solve_generic_dynamical_U.out.h5",
           "solve_generic_dynamical_U.ref.h5", precision=2e-9)
