# Test for a simple two-orbital system with no interaction. 
# We consider a hybridization function specified as either 
# 4 blocks of size 1 or 2 blocks of size 2. Results should be identical (move_move is disabled). 

from triqs.gf import *
from triqs.gf.descriptors import Function
from triqs.utility import mpi
from triqs.operators import n
from h5 import *
import numpy as np
from triqs.utility.h5diff import h5diff
from triqs_ctseg import SolverCore as Solver

# Numerical values 
n_orb = 2 # number of orbitals 
beta = 20.0 # inverse temperature
mu   = 0.0 # chemical potential
eps = 0.3 # hybridization levels
V = 0.5 # hybridization strengths
n_tau = 1001

# Solve parameters
solve_params = {
    "h_int": 0*n("up", 0)*n("down", 0),
    "chemical_potential": [mu] * 2 * n_orb,
    "length_cycle": 50,
    "n_warmup_cycles": 1000,
    "n_cycles": 10000,
    "measure_F_tau": True,
    "measure_nn_tau": True,
    "measure_nn_static": True,
    "move_move_segment": False
    }

# --------- 4 blocks of size 1 -----------

# Input structure: 4 blocks of size 1 
gf_struct = [("up1", 1), ("dn1", 1), ("up2", 1), ("dn2", 1)]

# Construct solver
S = Solver(beta = beta,
           n_tau = n_tau,
           gf_struct= gf_struct
           )

# Input Delta(tau)
for name, block in S.Delta_tau:
    Delta_iw = GfImFreq(indices=range(np.size(block)), beta=beta, n_points=n_tau//2)
    Delta_iw << V**2 * inverse(iOmega_n - eps)
    block << Fourier(Delta_iw)

# Solve
S.solve(**solve_params)

# Save output
if mpi.is_master_node(): 
    with HDFArchive("two_orbitals_4b.out.h5", "a") as A: 
        A['nn_tau'] = S.results.nn_tau
        A['nn'] = S.results.nn_static
        A['densities'] = S.results.densities

# --------- 2 blocks of size 2 -----------

# Input structure: 2 blocks of size 2 
gf_struct = [("up1", 2), ("dn1", 2)]

# Construct solver
S = Solver(beta = beta,
           n_tau = n_tau,
           gf_struct= gf_struct
           )

# Input Delta(tau)
delta_iw = GfImFreq(indices=[0, 1], beta=beta, n_points=n_tau//2)
delta_iw[0, 0] << V**2 * inverse(iOmega_n - eps)
delta_iw[1, 1] << V**2 * inverse(iOmega_n - eps)
Delta_iw = BlockGf(block_list = [delta_iw, delta_iw], name_list = ["up1", "dn1"])

S.Delta_tau << Fourier(Delta_iw)

# Solve
S.solve(**solve_params)

# Save output
if mpi.is_master_node(): 
    with HDFArchive("two_orbitals_2b.out.h5", "a") as A: 
        A['nn_tau'] = S.results.nn_tau
        A['nn'] = S.results.nn_static
        A['densities'] = S.results.densities

# --------- Compare outputs -----------
if mpi.is_master_node():
    h5diff("two_orbitals_2b.out.h5", "two_orbitals_4b.out.h5", precision=1e-9)