# Single orbital with Jperp interaction and no hybridization 
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
eps = 0.2
n_tau = 2051
n_tau_bosonic = 2001
wp = 1
g = 4

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

# Jperp interaction
d0_iw = GfImFreq(indices=[0], beta=beta, statistic="Boson", n_points=n_tau_bosonic//2)
d0_tau = GfImTime(indices=[0], beta=beta, statistic="Boson", n_points=n_tau_bosonic)
d0_iw << Function(lambda w: g * wp**2 / (w**2 - wp**2))
d0_tau << Fourier(d0_iw)
S.Jperp_tau << d0_tau

# Solve parameters
solve_params = {
    "h_int": U*n("up", 0)*n("down", 0),
    "h_loc0": -mu * (n("up", 0) + n("down", 0)),
    "length_cycle": 50,
    "n_warmup_cycles": 1000,
    "n_cycles": 10000,
    "measure_nn_tau": True,
    "measure_nn_static": True
    }

# Solve
S.solve(**solve_params)

# Save and compare to reference
if mpi.is_master_node():
    with h5.HDFArchive("Jperp.out.h5", 'w') as A:
        A['G_tau'] = S.results.G_tau
        A['nn_tau'] = S.results.nn_tau
        A['nn'] = S.results.nn_static
        A['densities'] = S.results.densities

    h5diff("Jperp.out.h5", "Jperp.ref.h5", precision=1e-9) 
