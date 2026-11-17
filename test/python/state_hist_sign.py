# Copyright (c) 2026, The Simons Foundation
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

"""The state histogram is weighted with the Monte-Carlo sign.

A hybridization Delta(tau) > 0 (not causal) gives configurations of both signs.
The histogram must still be normalized to one, and the occupations it implies
must agree with the separately measured, sign-weighted densities.
"""

import numpy as np
import triqs.utility.mpi as mpi
from triqs.operators import n
from triqs_ctseg import Solver

beta, U, mu = 10., 2., 1.
gf_struct = [('down', 1), ('up', 1)]

S = Solver(gf_struct=gf_struct, beta=beta, n_tau=201)
for _, d in S.Delta_tau:
    d.data[:] = 0.3

S.solve(h_int=U * n('up', 0) * n('down', 0), h_loc0=-mu * (n('up', 0) + n('down', 0)),
        length_cycle=20, n_warmup_cycles=500, n_cycles=10000, random_seed=23,
        measure_densities=True, measure_state_hist=True, measure_average_sign=True, measure_G_tau=False)

if mpi.is_master_node():
    sign = S.results.average_sign
    hist = np.asarray(S.results.state_hist)
    occ_hist = [hist[1] + hist[3], hist[2] + hist[3]]   # colors: down = bit 0, up = bit 1
    occ_meas = [S.results.densities['down'][0], S.results.densities['up'][0]]
    print(f"average sign {sign:.4f}, sum(state_hist) {hist.sum():.12f}")
    print(f"occupations from histogram {occ_hist}, measured {occ_meas}")
    assert abs(sign) < 0.99, "the test needs a sign problem"
    assert abs(hist.sum() - 1) < 1e-10
    assert np.allclose(occ_hist, occ_meas, atol=1e-10)
