from triqs.gf import *
from h5 import *
from triqs.operators import *
from numpy import matrix
import triqs.utility.mpi as mpi
from triqs_ctseg import SolverCore as Solver

U = 2.0
beta = 10.0
mu = 1.0
eps0 = 0.0
l=1.0
w0=1.0

gf_struct = [("up",1), ("down",1)]
S = Solver(beta = beta,
           gf_struct = gf_struct,
           n_tau = 10000,
           n_tau_nn = 1000,
           )

for b1 in dict(gf_struct).keys():
  for b2 in dict(gf_struct).keys():
      S.D0_iw[b1+"|"+b2] << l**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))
S.G0_iw << inverse(iOmega_n + mu - inverse(iOmega_n - eps0))

h_int = U * n('up',0)*n('down',0)

S.solve(h_int=h_int,
        n_cycles  = 5000,
        length_cycle = 100,
        n_warmup_cycles = 1000,
        n_w_b_vertex = 2,
        n_w_f_vertex = 20,
        measure_nn=True,
        measure_nnw=True,
        measure_gw=True,
        measure_fw=True,
        measure_g2w=True,
        measure_f2w=True,
        measure_g3w=True,
        measure_f3w=True,
        evaluate_vertex = True
        )

if mpi.is_master_node():
  A = HDFArchive("ct_seg_single_bath_vertex.out.h5",'a')
  A['G_tau'] = S.G_tau
  A['nn_iw'] = S.nn_iw
  A['G_iw']  = S.G_iw
  A['F_iw']  = S.F_iw
  A["G_2w"] = S.G_2w
  A["F_2w"] = S.F_2w
  A["G_3w"] = S.G_3w
  A["F_3w"] = S.F_3w
  del A
from triqs.utility.h5diff import h5diff
h5diff("ct_seg_single_bath_vertex.out.h5", "ct_seg_single_bath_vertex.ref.h5")

