import itertools, time
import numpy as np
from triqs.gf import *
from triqs.gf.descriptors import Function
from h5 import *
from triqs.operators import *
import triqs.utility.mpi as mpi

beta = 20.0
U = 2.0
mu = 0.2
eps0 = 0.3
V = 1.0

from triqs_ctseg import Solver
gf_struct = [("up",1), ("down",1)]
S = Solver(beta = beta,  gf_struct = gf_struct)
S.G0_iw << inverse(iOmega_n + mu - V**2 * inverse(iOmega_n - eps0))

w_ch, w_sp = 1.0, 0.8 #screening frequencies
sigz = lambda s1, s2 : 1 if s1==s2 else -1
vals, zeros = np.linspace(0,1,10), [0.0]*10

for g_ch, g_sp in list(zip(vals,zeros)) + list(zip(zeros,vals))[1:]:
 print(g_ch, g_sp)

 for b1, b2 in itertools.product(list(dict(gf_struct).keys()), list(dict(gf_struct).keys())):
  S.D0_iw[b1+"|"+b2] << Function(lambda w : g_ch**2 *2*w_ch/(w**2-w_ch**2) \
                                          + g_sp**2 *2*w_sp/(w**2-w_sp**2)*sigz(b1,b2) )
 S.Jperp_iw << Function(lambda w : 4 * g_sp**2 * 2*w_sp/(w**2-w_sp**2))

 L = [] 
 for i in range(5):
  t0 = time.clock()
  S.solve(h_int = U * n('up',0)*n('down',0),
         n_cycles = 50000, length_cycle = 100, n_warmup_cycles = 100)
  t1 = time.clock()
  L.append(t1 -t0)
  print("time =", t1 - t0)

 if mpi.is_master_node():
  with HDFArchive("dyn_interactions_g_ch_%s_g_sp_%s.h5"%(g_ch,g_sp),'w') as A: 
   A['G_tau'] = S.G_tau
   A["time"] = np.mean(L)
