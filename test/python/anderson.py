# Single orbital with static interaction at half-filling. 
from triqs.gf import *
import triqs.utility.mpi as mpi
from triqs.gf.descriptors import Function
from triqs.operators import n
from triqs.gf import iOmega_n
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

# Hybridization Delta(tau)
iw_mesh = MeshImFreq(beta, 'Fermion', n_tau//2)
delta_iw = GfImFreq(indices = [0], mesh = iw_mesh)
delta_iw << inverse(iOmega_n - eps)
S.Delta_tau["up"] << Fourier(delta_iw)
S.Delta_tau["down"] << Fourier(delta_iw)

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

# Save and compare to reference
if mpi.is_master_node():
    with h5.HDFArchive("anderson.out.h5", 'w') as A:
        A['G_tau'] = S.results.G_tau
        A['F_tau'] = S.results.F_tau
        A['nn_tau'] = S.results.nn_tau
        A['nn'] = S.results.nn_static
        A['densities'] = S.results.densities

    h5diff("anderson.out.h5", "anderson.ref.h5", precision=1e-9) 
