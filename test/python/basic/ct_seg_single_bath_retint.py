from triqs.gf import *
from triqs.gf.descriptors import Function
from h5 import *
from triqs.operators import *
from numpy import array,zeros, sinh,cosh, matrix
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

S = Solver(beta = beta,
           gf_struct = gf_struct,
           n_tau = 10000,
           )

h_int = U * n('up',0)*n('down',0)

l=1.0  #electron boson coupling
w0=1.0 #screening frequency

for b1 in dict(gf_struct).keys():
  for b2 in dict(gf_struct).keys():
      S.D0_iw[b1+"|"+b2] << Function(lambda w: l**2*2*w0/(w**2-w0**2))

# init the Green function
S.G0_iw << inverse(iOmega_n + mu - V**2 * inverse(iOmega_n - eps0))

S.solve(h_int=h_int,
        n_cycles  = 1000000,
        length_cycle = 1,
        n_warmup_cycles = 0,
        measure_nn=True)

if mpi.is_master_node():
  A = HDFArchive("ct_seg_single_bath_retint.out.h5",'w')
  A['G_tau'] = S.G_tau
  A['K_tau'] = S.K_tau
  A['Kprime_tau'] = S.Kprime_tau
from triqs.utility.h5diff import h5diff
h5diff("ct_seg_single_bath_retint.out.h5", "ct_seg_single_bath_retint.ref.h5")
