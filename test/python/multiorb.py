# Multi-orbital impurity with diagonal Delta(tau). 
# Number of orbitals can be changed. 
from triqs.gf import *
from triqs.gf.descriptors import Function
from triqs.utility import mpi
from triqs.operators import *
from h5 import *
import numpy as np
from triqs_ctseg import SolverCore as Solver
import h5
from triqs.utility.h5diff import h5diff

# -------------- Model parameters -----------

beta = 20.0 # inverse temperature
mu   = 0.2 # chemical potential

#number of orbitals
n_orb = 3

#interaction
U    = 2.0
Up   = 1.0
J    = 0.0

#Hamiltonian for the model
h_int = 0
for l in range(n_orb):
    h_int += U * n(f'up{l+1}', 0) * n(f'down{l+1}', 0)
for l1 in range(n_orb):
    for l2 in range(l1+1, n_orb):
        h_int += (Up - J) * (n(f'up{l1+1}', 0) * n(f'up{l2+1}', 0) + n(f'down{l1+1}', 0) * n(f'down{l2+1}', 0))
        h_int += Up * (n(f'up{l1+1}', 0) * n(f'down{l2+1}', 0) + n(f'down{l1+1}', 0) * n(f'up{l2+1}', 0))

#hybridization levels
eps ={} 
for l in range(n_orb):
    eps[f'up{l+1}'] = 0.3
    eps[f'down{l+1}'] = 0.3

#hybridization strengths
V ={}
for l in range(n_orb):
    V[f'up{l+1}'] = 1./np.sqrt(l+1)
    V[f'down{l+1}'] = 1./np.sqrt(l+1)


# -------------- Construct solver -----------
gf_struct = [(f'up{(l+2)//2}', 1) if l % 2 == 0 else (f'down{(l+1)//2}', 1) for l in range(2 * n_orb)]
S = Solver(beta = beta,
           n_tau = 10000,
           gf_struct= gf_struct
           )

# -------------- Input Delta(tau) -----------
for name, block in S.Delta_tau:
    Delta_iw = GfImFreq(indices=[0], beta=beta, n_points=1000)
    Delta_iw << V[name]**2 * inverse(iOmega_n - eps[name])
    block << Fourier(Delta_iw)

# -------------- Solve -----------
S.solve(h_int=h_int, 
        hartree_shift=[mu] * 2 * n_orb,
        n_cycles  = 10000,
        length_cycle = 100,
        n_warmup_cycles = 1000,     
        measure_gt=True,
        measure_nn=False,
        )

if mpi.is_master_node():
    with h5.HDFArchive("multiorb.out.h5", 'w') as A:
        A['G_tau'] = S.results.G_tau

# --------- Compare to reference ----------      
    h5diff("multiorb.out.h5", "multiorb.ref.h5", precision=1e-9)