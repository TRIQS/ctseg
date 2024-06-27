# Single orbital with dynamical spin-spin interactions. 
# Can be compared with reference obtained with CTINT (ctint.ref.h5)
from triqs.gf import *
import triqs.utility.mpi as mpi
from triqs.gf.descriptors import Function
from triqs.operators import n
import h5
from triqs.utility.h5diff import h5diff
from triqs_ctseg import SolverCore as Solver

# Numerical values
beta = 10
U = 4.0
mu = U/2
J = 0.5
n_tau = 2051
n_tau_bosonic = 2001

# Solver construction parameters
gf_struct = [('down', 1), ('up', 1)]
constr_params = {
    "gf_struct": gf_struct,
    "beta": beta,
    "n_tau": n_tau,
    "n_tau_bosonic": n_tau_bosonic
}

# Construct solver
S = Solver(**constr_params)

# Get inputs from reference file
with h5.HDFArchive("ctint.ref.h5", 'r') as Af:
    g0 = Af["dmft_loop/i_001/S/G0_iw/up"]
    Q_tau = Af["dmft_loop/i_000/Q_tau"]

# Hybridization Delta(tau)
n_iw = 1025
Delta = GfImFreq(indices=[0], beta=beta, n_points=n_iw)
invg0 = GfImFreq(indices=[0], beta=beta, n_points=n_iw)
invg0 << inverse(g0)
Delta << iOmega_n + mu - invg0
S.Delta_tau << Fourier(Delta)

# Spin-spin interaction (D0(tau) and Jperp(tau))
S.Jperp_tau << -J**2*Q_tau
S.D0_tau["up", "up"] << -0.25*J**2*Q_tau
S.D0_tau["down", "down"] << -0.25*J**2*Q_tau
S.D0_tau["up", "down"] << 0.25*J**2*Q_tau
S.D0_tau["down", "up"] << 0.25*J**2*Q_tau

# Solve parameters
solve_params = {
    "h_int": U*n("up", 0)*n("down", 0),
    "h_loc0": -mu * (n("up", 0) + n("down", 0)),
    "length_cycle": 50,
    "n_warmup_cycles": 1000,
    "n_cycles": 10000,
    "measure_F_tau": True,
    "measure_nn_tau": True,
    "measure_nn_static": True
    }

# Solve
S.solve(**solve_params)

# Unwrap Block2Gf nn_tau into matrix Gf
n_color = 2
nn_tau = GfImTime(indices=[0, 1], beta=beta, statistic = "Boson", n_points=n_tau_bosonic)
nn_tau[0, 0] = S.results.nn_tau["down", "down"][0, 0]
nn_tau[1, 0] = S.results.nn_tau["up", "down"][0, 0]
nn_tau[0, 1] = S.results.nn_tau["down", "up"][0, 0]
nn_tau[1, 1] = S.results.nn_tau["up", "up"][0, 0]

# Save and compare to reference
if mpi.is_master_node():
    with h5.HDFArchive("spin_spin.out.h5", 'w') as A:
        A['G_tau'] = S.results.G_tau
        A['F_tau'] = S.results.F_tau
        A['nn_tau'] = nn_tau
        A['nn'] = S.results.nn_static
        A['densities'] = S.results.densities

    h5diff("spin_spin.out.h5", "spin_spin.ref.h5", precision=1e-9) 
