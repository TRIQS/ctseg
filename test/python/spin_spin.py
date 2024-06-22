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
n_tau_k = 2001

# Solver construction parameters
constr_params = {
    "gf_struct": [('down', 1), ('up', 1)],
    "beta": beta,
    "n_tau": n_tau,
    "n_tau_k": n_tau_k
}

# Construct solver
S = Solver(**constr_params)

# Get inputs from reference file
with h5.HDFArchive("ctint.ref.h5", 'r') as Af:
    g0 = Af["dmft_loop/i_001/S/G0_iw/up"]
    Q_tau = Af["dmft_loop/i_000/Q_tau"]

# Hybridization Delta(tau)
n_iw = 1025
delta = GfImFreq(indices=[0], beta=beta, n_points=n_iw)
invg0 = GfImFreq(indices=[0], beta=beta, n_points=n_iw)
invg0 << inverse(g0)
delta << iOmega_n + mu - invg0
S.Delta_tau << Fourier(delta)

# Spin-spin interaction (D0(tau) and Jperp(tau))
D0 = GfImTime(indices=[0, 1], beta=beta, statistic='Boson', n_points=n_tau_k)
D0[0, 0] = -Q_tau[0, 0]
D0[1, 1] = -Q_tau[0, 0]
D0[0, 1] = Q_tau[0, 0]
D0[1, 0] = Q_tau[0, 0]
S.Jperp_tau << -J**2*Q_tau
S.D0_tau << 0.25*J**2*D0

# Solve parameters
solve_params = {
    "h_int": U*n("up", 0)*n("down", 0),
    "chemical_potential": [mu, mu],
    "length_cycle": 50,
    "n_warmup_cycles": 1000,
    "n_cycles": 10000,
    "measure_F_tau": True,
    "measure_nn_tau": True,
    "measure_nn_static": True
    }

# Solve
S.solve(**solve_params)

# Save and compare to reference
if mpi.is_master_node():
    with h5.HDFArchive("spin_spin.out.h5", 'w') as A:
        A['G_tau'] = S.results.G_tau
        A['F_tau'] = S.results.F_tau
        A['nn_tau'] = S.results.nn_tau
        A['nn'] = S.results.nn_static
        A['densities'] = S.results.densities

    h5diff("spin_spin.out.h5", "spin_spin.ref.h5", precision=1e-9) 
