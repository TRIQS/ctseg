# Copyright (c) 2026, The Simons Foundation
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

"""End-to-end test of solve_generic on MeshImFreq with static interaction.

Mirrors the physics of anderson.py so the tau-space observables are
physically comparable.  Exercises:
  - _prepare_delta_tau MeshImFreq branch (tail-fit Fourier)
  - postprocess_sigma improved-estimator path (F_tau is always
    measured by solve_generic)
  - tail-fit Sigma_HF path (analytic_hf=False, the default)
  - h_loc0 supplied as a list of block matrices (one of two accepted forms)
"""

import numpy as np
from triqs.gfs import MeshImFreq, Gf, BlockGf, inverse, iOmega_n
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
n_iw = 1025
n_tau = 2051

# Hybridization on MeshImFreq
iw_mesh = MeshImFreq(beta=beta, S="Fermion", n_iw=n_iw)
gf_struct = [('down', 1), ('up', 1)]
Delta_iw = BlockGf(
    name_block_generator=[(bl, Gf(mesh=iw_mesh, target_shape=[dim, dim]))
                          for bl, dim in gf_struct],
    make_copies=False,
)
Delta_iw << inverse(iOmega_n - eps)

# h_loc0 passed as a list of dense block matrices (one per block in gf_struct).
# The Operator form is exercised by the other solve_generic_* tests.
h_loc0 = [-mu * np.eye(dim) for _, dim in gf_struct]
h_int = U * n('up', 0) * n('down', 0)

results = solve_generic(
    Delta_iw, h_loc0, h_int,
    n_iw=n_iw,
    n_tau=n_tau,
    length_cycle=50,
    n_warmup_cycles=1000,
    n_cycles=10000,
)

if mpi.is_master_node():
    # postprocess_pi does not run (solver.results.nn_nu is not measured by
    # ctseg), so Pi_iw/W_iw/Chi_iw should not appear in results
    assert 'Pi_iw' not in results.keys()
    assert results.G_iw is not None
    assert results.Sigma_iw is not None
    assert results.Sigma_dynamic is not None
    assert results.Sigma_HartreeFock is not None
    assert results.G_tau is not None

    with h5.HDFArchive("solve_generic_imfreq.out.h5", 'w') as A:
        A['G_tau'] = results.G_tau
        A['G_iw'] = results.G_iw
        A['Sigma_iw'] = results.Sigma_iw
        A['Sigma_dynamic'] = results.Sigma_dynamic
        A['Sigma_HartreeFock'] = results.Sigma_HartreeFock
        A['densities'] = results.Solver.results.densities

    h5diff("solve_generic_imfreq.out.h5", "solve_generic_imfreq.ref.h5", precision=1e-9)
