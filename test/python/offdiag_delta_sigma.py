# Copyright (c) 2026, The Simons Foundation
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

# Two orbitals per spin coupled to four bath sites, so that Delta(iw), G(iw) and
# Sigma(iw) have sizeable off-diagonal elements within each block.
# The improved-estimator self-energy returned by solve_generic (Sigma = F G^{-1})
# must agree with the Dyson self-energy at low Matsubara frequencies, where the
# Dyson equation is still accurate.

import numpy as np

import triqs.utility.mpi as mpi
from triqs.gfs import BlockGf, MeshImFreq, inverse, make_hermitian
from triqs.gfs.tools import make_zero_tail
from triqs.operators.util import U_matrix_kanamori, h_int_density

from triqs_ctseg.solve_generic import solve_generic

beta = 10.0
n_orb = 2
U, J = 2.0, 0.3
mu = 2.0
eps_imp = np.array([0.0, 0.2])
eps_bath = np.array([-1.2, -0.4, 0.4, 1.2])
V = np.array([[0.6, 0.6, 0.6, 0.6],
              [0.6, 0.2, 0.6, 0.2]])
n_iw = 200
n_tau = 2001
n_check = 2          # number of positive Matsubara frequencies compared
tol = 0.2            # QMC noise of the difference is < 0.1; the element-wise F/G gives ~1.5

block_names = ['up', 'dn']
gf_struct = [(bl, n_orb) for bl in block_names]

Umat, Upmat = U_matrix_kanamori(n_orb, U_int=U, J_hund=J)
h_int = h_int_density(block_names, n_orb, Umat, Upmat, off_diag=True)
h_loc0_mat = np.diag(eps_imp - mu)

mesh = MeshImFreq(beta=beta, statistic='Fermion', n_iw=n_iw)
Delta_iw = BlockGf(mesh=mesh, gf_struct=gf_struct)
G0_iw = Delta_iw.copy()
for bl in block_names:
    for iw in mesh:
        delta = (V / (iw.value - eps_bath)) @ V.T
        Delta_iw[bl][iw] = delta
        G0_iw[bl][iw] = np.linalg.inv(iw.value * np.eye(n_orb) - h_loc0_mat - delta)

res = solve_generic(
    Delta_iw, [h_loc0_mat] * len(block_names), h_int,
    n_iw=n_iw, n_tau=n_tau,
    n_warmup_cycles=2000, n_cycles=40000, length_cycle=50,
    random_seed=34788 + 928374 * mpi.rank,
)

G_iw = Delta_iw.copy()
for bl in block_names:
    known_moments = make_zero_tail(G_iw[bl], n_moments=2)
    known_moments[1] = np.eye(n_orb)
    G_iw[bl].set_from_fourier(res.Solver.results.G_tau[bl], known_moments)
G_iw << make_hermitian(G_iw)
Sigma_dyson = inverse(G0_iw) - inverse(G_iw)

if mpi.is_master_node():
    # the model is spin symmetric: average the blocks to reduce the noise
    s_ie = sum(res.Sigma_iw[bl].data[n_iw:n_iw + n_check] for bl in block_names) / len(block_names)
    s_dy = sum(Sigma_dyson[bl].data[n_iw:n_iw + n_check] for bl in block_names) / len(block_names)
    # the off-diagonal self-energy must be non-trivial for the test to be meaningful
    assert np.abs(s_dy[0, 0, 1]) > 0.1, s_dy[0]
    diff = np.abs(s_ie - s_dy).max()
    print(f"max |Sigma_IE - Sigma_Dyson| over {n_check} frequencies = {diff:.4f}")
    assert diff < tol, f"improved estimator and Dyson self-energies differ by {diff}"
