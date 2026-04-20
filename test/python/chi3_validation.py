# Copyright (c) 2024--present, The Simons Foundation
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

# Validate chi3 (DLR2D, NFFT) against g2w (uniform mesh, naive FT).
# Both use the same nw[0] = beta*(n-1) convention, so they agree at all
# bosonic frequencies including Omega = 0.

from triqs.gfs import *
from triqs.operators import n
import triqs.utility.mpi as mpi
import numpy as np
from triqs_ctseg import SolverCore as Solver

# Parameters
beta = 10
U = 4.0
mu = U / 2
eps = 0.2
n_tau = 2051
n_tau_bosonic = 2001

# Solver construction
gf_struct = [('down', 1), ('up', 1)]
S = Solver(gf_struct=gf_struct, beta=beta, n_tau=n_tau, n_tau_bosonic=n_tau_bosonic)

# Hybridization
iw_mesh = MeshImFreq(beta, 'Fermion', n_tau // 2)
delta_iw = GfImFreq(indices=[0], mesh=iw_mesh)
delta_iw << inverse(iOmega_n - eps)
S.Delta_tau["up"] << Fourier(delta_iw)
S.Delta_tau["down"] << Fourier(delta_iw)

# Solve with both g2w and chi3
n_w_b = 4
n_w_f = 5
solve_params = {
    "h_int": U * n("up", 0) * n("down", 0),
    "h_loc0": -mu * (n("up", 0) + n("down", 0)),
    "length_cycle": 50,
    "n_warmup_cycles": 2000,
    "n_cycles": 40000,
    "random_seed": 12345,
    # g2w (reference)
    "measure_g2w": True,
    "n_w_b_vertex": n_w_b,
    "n_w_f_vertex": n_w_f,
    # chi3 via NFFT on DLR2D mesh
    "measure_chi3": True,
    "dlr_wmax": 10.0,
    "dlr_eps": 1e-10,
}
S.solve(**solve_params)

if mpi.is_master_node():
    g2w = S.results.g2w
    chi3 = S.results.chi3

    # Convert chi3 from DLR2D to DLR2D coefficient form (evaluable at arbitrary frequencies)
    chi3_dlr2d = {bl: make_gf_dlr2d(chi3[bl]) for bl in chi3.indices}

    # Compare: for each (Omega, nu) point in the g2w mesh, evaluate chi3 at the
    # corresponding fermionic pair (nu1, nu2) = (nu, Omega + nu) [PH channel]
    mesh_b = MeshImFreq(beta, 'Boson', n_w_b)
    mesh_f = MeshImFreq(beta, 'Fermion', n_w_f)

    max_diff = 0.0
    for bl in g2w.indices:
        b1, b2 = bl
        bl1_size = g2w[bl].target_shape[0]
        bl2_size = g2w[bl].target_shape[2]
        for iW in mesh_b:
            for inu in mesh_f:
                Omega = iW.value
                nu = inu.value
                ferm_ph = (nu, Omega + nu)  # PH channel: (nu1, nu2) = (nu, Omega + nu)
                g2w_val = g2w[bl][iW, inu]
                chi3_val = chi3_dlr2d[bl](ferm_ph)
                diff = np.max(np.abs(g2w_val - chi3_val))
                max_diff = max(max_diff, diff)

    print(f"Max |g2w - chi3| = {max_diff:.2e}")
    assert max_diff < 1e-4, f"chi3 vs g2w mismatch: max diff = {max_diff:.2e}"
    print("PASSED: chi3 agrees with g2w")
