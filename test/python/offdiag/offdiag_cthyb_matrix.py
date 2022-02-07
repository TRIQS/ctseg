from triqs.gf import *
from triqs.operators import *
from cthyb_matrix import SolverCore as Solver

# Parameters
V = 0.5
U = 3.0
mu = U/2.
beta = 10
e_h=0.5
alpha=0.5
split=.05

# Construct the impurity solver with the inverse temperature
# and the structure of the Green's functions
S = Solver(beta = beta, gf_struct = [ ('up',2) ])

# Initialize the non-interacting Green's function S.G0
for spin, g0 in S.G0 :
    g0[0,0] = iOmega_n + mu + split/2 - V**2 * (inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h)) 
    g0[1,0] =  - alpha*V**2 *(inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h))   
    g0[1,1] = iOmega_n + mu - split/2 - V**2 *(inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h)) 
    g0[0,1] = - alpha* V**2 *(inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h))
    g0.invert()
# Run the solver. The result will be in S.G
S.solve(H_local = U * N('up',1) * N('up',2),   # Local Hamiltonian 
        quantum_numbers = {                      # Quantum Numbers 
          'Nup1' : N('up',1),                     # Operators commuting with H_Local
          'Nup2' : N('up',2) },          
        n_cycles  = 500000,                      # Number of QMC cycles
        length_cycle = 100,                      # Length of one cycle 
        n_warmup_cycles = 1000,                 # Warmup cycles
        n_legendre = 50,                         # Number of Legendre coefficients
        random_name = 'mt19937',                 # Name of the random number generator
        use_segment_picture = False,              # Use the segment picture
        measured_operators = {                   # Operators to be averaged
          'Nimp' : N('up',1)+N('up',2) }
        )

# Save the results in an hdf5 file (only on the master node)
from h5 import HDFArchive
import triqs.utility.mpi as mpi

if mpi.is_master_node():
  Results = HDFArchive("example_2.h5",'w')
  Results["G"] = S.G
  Results["Gl"] = S.G_legendre
  Results["Nimp"] = S.measured_operators_results['Nimp']
  Results["G0"] = S.G0
  Results["Delta_tau"] = S.Delta_tau
