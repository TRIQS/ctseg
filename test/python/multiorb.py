# Multi-orbital impurity with diagonal Delta(tau). 
# Number of orbitals can be changed. 
from triqs.gf import *
from triqs.gf.descriptors import Function
from triqs.utility import mpi
from triqs.operators import n
from h5 import *
import numpy as np
from triqs.utility.h5diff import h5diff
from triqs_ctseg import SolverCore as Solver

# Number of orbitals
n_orb = 3

# Numerical values
beta = 20.0 # inverse temperature
mu   = 0.2 # chemical potential
eps = 0.3 # hybridization levels
V = 1. # hybridization strengths
U    = 2.0 # same orbital interaction 
Up   = 1.0 # different orbital interaction 
J    = 0.1 # Hund coupling
n_tau = 1001
n_tau_k = 1001

# Solver construction parameters
constr_params = {
    "gf_struct": [(f'up{(l+2)//2}', 1) if l % 2 == 0 else (f'down{(l+1)//2}', 1) for l in range(2 * n_orb)],
    "beta": beta,
    "n_tau": n_tau,
    "n_tau_k": n_tau_k
}

# Construct solver
S = Solver(**constr_params)

# Interaction Hamiltonian
h_int = 0
for l in range(n_orb):
    h_int += U * n(f'up{l+1}', 0) * n(f'down{l+1}', 0)
for l1 in range(n_orb):
    for l2 in range(l1+1, n_orb):
        h_int += (Up - J) * (n(f'up{l1+1}', 0) * n(f'up{l2+1}', 0) + n(f'down{l1+1}', 0) * n(f'down{l2+1}', 0))
        h_int += Up * (n(f'up{l1+1}', 0) * n(f'down{l2+1}', 0) + n(f'down{l1+1}', 0) * n(f'up{l2+1}', 0))

# Hybridization Delta(tau)
for name, block in S.Delta_tau:
    Delta_iw = GfImFreq(indices=[0], beta=beta, n_points=n_tau//2)
    Delta_iw << V**2 * inverse(iOmega_n - eps)
    block << Fourier(Delta_iw)

# Solve parameters
solve_params = {
    "h_int": h_int,
    "hartree_shift": [mu] * 2 * n_orb,
    "length_cycle": 50,
    "n_warmup_cycles": 1000,
    "n_cycles": 10000,
    "measure_ft": True,
    "measure_nnt": True,
    "measure_nn": True
    }

# Solve
S.solve(**solve_params)

# Save and compare to reference
if mpi.is_master_node():
    with HDFArchive("multiorb.out.h5", 'w') as A:
        A['G_tau'] = S.results.G_tau
        A['F_tau'] = S.results.F_tau
        A['nn_tau'] = S.results.nn_tau
        A['nn'] = S.results.nn_static
        A['densities'] = S.results.densities

    #h5diff("multiorb.out.h5", "multiorb.ref.h5", precision=1e-9)