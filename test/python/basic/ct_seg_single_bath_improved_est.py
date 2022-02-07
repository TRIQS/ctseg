from triqs.gf import *
from h5 import *
from numpy import matrix
from triqs.operators import *
import triqs.utility.mpi as mpi

# set up a few parameters
U = 2.0
Beta = 20.0
mu = 0.2
eps0 = 0.3
V = 1.0

# Construct a Segment solver
from triqs_ctseg import SolverCore as Solver

S = Solver(beta = Beta,                                                      # inverse temperature
           gf_struct = [("up",1),("down",1)],
           n_tau = 10000,
           )

h_int = U * n('up',0)*n('down',0)

# init the Green function
S.G0_iw << inverse(iOmega_n + mu - V**2 * inverse(iOmega_n - eps0))

S.solve(h_int=h_int,
        measure_gt=False,
        measure_gw=True,
        measure_fw=True,
        n_cycles  = 3000,
        length_cycle = 50,
        n_warmup_cycles = 0,
        measure_nn=True)

if mpi.is_master_node():
  A = HDFArchive("ct_seg_single_bath_improved_est.out.h5",'w')
  A['G_iw'] = S.G_iw
  A['F_iw'] = S.F_iw
  A['Sigma_iw'] = S.Sigma_iw
from triqs.utility.h5diff import h5diff
h5diff("ct_seg_single_bath_improved_est.out.h5", "ct_seg_single_bath_improved_est.ref.h5")
