# Copyright (c) 2026, The Simons Foundation
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

"""End-to-end test of solve_density on MeshImFreq with static interaction.

solve_density is the density-only variant of solve_generic: it turns
off measure_G_tau / measure_F_tau / measure_nn_tau and skips the
postprocess pipeline.  The only output of interest is the converged
block-resolved density vector.
"""

from triqs.gf import MeshImFreq, Gf, BlockGf, inverse, iOmega_n
from triqs.operators import n
import triqs.utility.mpi as mpi
import h5
from triqs.utility.h5diff import h5diff

from triqs_ctseg import solve_density

# Model
beta = 10.
U = 4.
mu = U / 2.
eps = 0.2
n_iw = 1025
n_tau = 2051

iw_mesh = MeshImFreq(beta=beta, S='Fermion', n_iw=n_iw)
gf_struct = [('down', 1), ('up', 1)]
Delta_iw = BlockGf(
    name_block_generator=[(bl, Gf(mesh=iw_mesh, target_shape=[dim, dim]))
                          for bl, dim in gf_struct],
    make_copies=False,
)
Delta_iw << inverse(iOmega_n - eps)

h_loc0 = -mu * (n('up', 0) + n('down', 0))
h_int = U * n('up', 0) * n('down', 0)

results = solve_density(
    Delta_iw, h_loc0, h_int,
    n_iw=n_iw,
    n_tau=n_tau,
    length_cycle=50,
    n_warmup_cycles=1000,
    n_cycles=10000,
)

if mpi.is_master_node():
    # density-only mode: only densities are populated, no postprocessing ran
    assert results.Solver is not None
    assert results.Solver.results.densities is not None
    assert results.G_iw is None
    assert results.Sigma_iw is None
    assert results.Sigma_HartreeFock is None

    with h5.HDFArchive("solve_density.out.h5", 'w') as A:
        A['densities'] = results.Solver.results.densities

    h5diff("solve_density.out.h5", "solve_density.ref.h5", precision=1e-9)
