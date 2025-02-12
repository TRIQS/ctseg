# Retarded interaction with move_double enabled. 
from triqs.gf import *
from triqs.gf.descriptors import Function
from triqs.utility import mpi
from triqs.operators import n
from h5 import *
import numpy as np
from triqs.utility.h5diff import h5diff
from triqs_ctseg import Solver

# Numerical values
beta = 10.0 # inverse temperature
mu   = -0.5 # chemical potential
wp   = 1.0 # Phonon frequency
n_tau = 2051
n_tau_bosonic = 2051

# Solver construction parameters
gf_struct = [("up", 1), ("down", 1)]
constr_params = {
    "gf_struct": gf_struct,
    "beta": beta,
    "n_tau": n_tau,
    "n_tau_bosonic": n_tau_bosonic
}

# Construct solver
S = Solver(**constr_params)

# Interaction Hamiltonian
h_int = 0 * n("up", 0) * n("down", 0)

# Local Hamiltonian
h_loc0 = -mu * (n("up", 0) + n("down", 0))

# Hybridization Delta(tau) and retarded D0(tau)
n_iw = 1025
Delta = GfImFreq(indices=[0], beta=beta, n_points=n_iw)
D0    = GfImFreq(indices=[0], beta=beta, n_points=n_iw, statistic="Boson")
Delta << inverse(iOmega_n + mu)
freq_boson = 1.j * np.array([(2*p)*np.pi/beta for p in range(-n_iw+1, n_iw)])
D0.data[:, 0, 0] = 1 / (freq_boson**2 - wp**2)
S.Delta_tau << Fourier(Delta)
S.D0_tau["up", "up"] << Fourier(D0)
S.D0_tau["down", "down"] << Fourier(D0)
S.D0_tau["up", "down"] << Fourier(D0)
S.D0_tau["down", "up"] << Fourier(D0)

# Solve parameters
solve_params = {
    "h_int": h_int,
    "h_loc0": h_loc0,
    "length_cycle": 50,
    "n_warmup_cycles": 1000,
    "n_cycles": 100000,
    "measure_F_tau": True,
    "measure_nn_tau": True,
    "measure_nn_static": True,
    "move_double_insert_segment": True,
    "move_double_remove_segment": True,
    }

# Solve
S.solve(**solve_params)

# Save and compare to reference
if mpi.is_master_node():
    with HDFArchive("move_double.out.h5", 'w') as A:
        A['G_tau'] = S.results.G_tau
        A['F_tau'] = S.results.F_tau
        A["nn_tau"] = S.results.nn_tau
        A['nn'] = S.results.nn_static
        A['densities'] = S.results.densities

    h5diff("move_double.out.h5", "move_double.ref.h5", precision=1e-9)
