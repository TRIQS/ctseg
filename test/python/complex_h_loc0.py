# Copyright (c) 2024--present, The Simons Foundation
# Copyright (c) 2024--present, Max Planck Institute for Polymer Research, Mainz, Germany
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

# h_loc0 built from complex matrices carries complex-flagged coefficients (real_or_complex)
# even when the imaginary part is numerical noise. CT-SEG uses a real local Hamiltonian:
#  - an imaginary part on the diagonal below imag_threshold must be dropped (taking the real
#    part), so a complex-typed h_loc0 reproduces the equivalent real one exactly;
#  - an imaginary part above imag_threshold must error.
from triqs.gfs import *
import triqs.utility.mpi as mpi
from triqs.operators import n
from triqs.gfs import iOmega_n
from triqs.utility.h5diff import h5diff
import h5
from triqs_ctseg import Solver

beta, U, mu, eps, n_tau, n_tau_bosonic = 10, 4.0, 2.0, 0.2, 2051, 2001
constr_params = {"gf_struct": [('down', 1), ('up', 1)], "beta": beta,
                 "n_tau": n_tau, "n_tau_bosonic": n_tau_bosonic}

def make_solver():
    S = Solver(**constr_params)
    iw_mesh = MeshImFreq(beta, 'Fermion', n_tau // 2)
    delta_iw = GfImFreq(indices=[0], mesh=iw_mesh)
    delta_iw << inverse(iOmega_n - eps)
    S.Delta_tau["up"] << Fourier(delta_iw)
    S.Delta_tau["down"] << Fourier(delta_iw)
    return S

base = {"h_int": U * n("up", 0) * n("down", 0), "length_cycle": 50,
        "n_warmup_cycles": 100, "n_cycles": 500, "random_seed": 42, "measure_F_tau": True}

# Real h_loc0
S_real = make_solver()
S_real.solve(h_loc0=-mu * (n("up", 0) + n("down", 0)), **base)

# Complex-typed h_loc0 with a negligible imaginary part: the imaginary part must be dropped,
# reproducing the real run (identical real coefficients, same seed -> same Markov chain).
S_cplx = make_solver()
S_cplx.solve(h_loc0=complex(-mu, 1e-15) * (n("up", 0) + n("down", 0)), **base)

if mpi.is_master_node():
    with h5.HDFArchive("complex_h_loc0_real.h5", 'w') as A:
        A['G_tau'] = S_real.results.G_tau
        A['F_tau'] = S_real.results.F_tau
        A['densities'] = S_real.results.densities
    with h5.HDFArchive("complex_h_loc0_cplx.h5", 'w') as A:
        A['G_tau'] = S_cplx.results.G_tau
        A['F_tau'] = S_cplx.results.F_tau
        A['densities'] = S_cplx.results.densities
    h5diff("complex_h_loc0_cplx.h5", "complex_h_loc0_real.h5", precision=1e-10)

# An imaginary part above imag_threshold (default 1e-13) must raise.
raised = False
try:
    S_bad = make_solver()
    S_bad.solve(h_loc0=complex(-mu, 1e-3) * (n("up", 0) + n("down", 0)), **base)
except RuntimeError:
    raised = True
assert raised, "an h_loc0 with Im = 1e-3 > imag_threshold should raise an error"

# imag_threshold is tunable: raising it above the imaginary part accepts it silently.
S_tuned = make_solver()
S_tuned.solve(h_loc0=complex(-mu, 1e-9) * (n("up", 0) + n("down", 0)), imag_threshold=1e-8, **base)
