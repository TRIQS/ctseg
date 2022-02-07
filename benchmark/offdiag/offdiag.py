from triqs.gf import *
from triqs.operators import *
from h5 import *
import triqs.utility.mpi as mpi
import time
import numpy as np

solver = "cthyb" 
solver = "ctseg" 

if solver == "ctseg":
 from triqs_ctseg import Solver
elif solver == "cthyb":
 from triqs.applications.impurity_solvers.cthyb import Solver

#for beta in [1., 2., 10., 20., 50., 100., 200., 500.]:
for beta in [100.]:
 U = 3.0
 mu = U/2.
 V = 0.5
 e_h = 0.5
 split = .5

 S = Solver(beta = beta, gf_struct = [("up",3), ("down",3)])

 alpha = {"up" : 0.25, "down": 0.5}
 E_bath = {sp: np.array([[-e_h, 0., alpha[sp]*e_h], [0., 0., 0.], [alpha[sp]*e_h, 0., e_h]]) for sp in alpha.keys()}
 E_loc = np.array([[split/2, 0., 0.],[0., 0., 0.],[0., 0., -split/2]]) 

 # initialize the Green function
 for spin in ['up','down']:
  S.G0_iw[spin] << iOmega_n + mu + E_loc - V**2 * inverse(iOmega_n - E_bath[spin])
 S.G0_iw.invert()

 L=[]
 for i in range(1):
  t0=time.clock()
  S.solve(h_int = U*sum([n(sp,i)*n(sp,j) for sp in ["up","down"] for i,j in [(0,1),(0,2),(1,2)]]),
         n_cycles = 500000, length_cycle = 500, n_warmup_cycles = 1000, move_move = False, move_swap_empty_lines = True )
  t1=time.clock()
  print("time = ", t1-t0)
  L.append(t1-t0)
 if mpi.is_master_node():
  with HDFArchive("offdiag_%s_beta_%s_three_bands.h5"%(solver,beta),'w') as A:
   A['G_tau'] = S.G_tau
   A['time'] = np.mean(L)
   print("avg time = ", np.mean(L))
