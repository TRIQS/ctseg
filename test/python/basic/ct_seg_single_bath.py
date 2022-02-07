from triqs.gf import *
from h5 import *
from numpy import matrix
from triqs.operators import *
import triqs.utility.mpi as mpi

# set up a few parameters
U = 2.0
beta = 20.0
mu = 0.2
eps0 = 0.3
V = 1.0

# Construct a Segment solver
from triqs_ctseg import SolverCore as Solver

gf_struct = [("up",1), ("down",1)]

S = Solver(beta = beta, # inverse temperature
           gf_struct = gf_struct,
           n_tau = 10000
           )

h_int = U * n('up',0)*n('down',0)

# init the Green function
S.G0_iw << inverse(iOmega_n + mu - V**2 * inverse(iOmega_n - eps0))

S.solve(
        h_int=h_int,
        move_move=False,
        move_swap_empty_lines=False,        
        n_cycles  = 100000,
        length_cycle = 1,
        n_warmup_cycles = 0,
        measure_nn=True,
        measure_nnt=True,
        measure_nnw=True,
        )

if mpi.is_master_node():
  A = HDFArchive("ct_seg_single_bath.out.h5",'w')
  A['G_tau'] = S.G_tau
  A['nn'] = S.nn
  A['nn_tau'] = S.nn_tau
  A['nn_iw'] = S.nn_iw
  A['G0'] = S.G0_iw
  print(S.average_sign)
  #A['mc_sign'] = S.average_sign
from triqs.utility.h5diff import h5diff
h5diff("ct_seg_single_bath.out.h5", "ct_seg_single_bath.ref.h5")
