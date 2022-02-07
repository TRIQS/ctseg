from triqs.gf import *
from triqs.gf.descriptors import Function
from triqs.operators import *
from h5 import *
import triqs.utility.mpi as mpi
from numpy import matrix,array,sinh,cosh  

# set up a few parameters
U = 1.  #Hubbard U
beta = 10.0  #inverse temperature
mu = .5 #chemical potential
V = .1 #hybridization strength
eps0=.2  #level of fermionic bath
g=0.3  #spin-boson coupling
w0=1.0 #bosonic energy

# Construct a Segment solver

gf_struct = [("up",1), ("down",1)]

from triqs_ctseg import SolverCore as Solver
S = Solver(beta = beta,  # inverse temperature
           gf_struct = gf_struct,
           n_tau = 10000,
           n_tau_nn = 10001
           )

h_int = U * n('up',0)*n('down',0)

sigz = {'up':{'up':1,'down':-1}, 'down':{'up':-1,'down':1}}

for b1 in dict(gf_struct).keys():
 for b2 in dict(gf_struct).keys():
   S.D0_iw[b1+"|"+b2] << g**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))*sigz[b1][b2]*0.25
S.Jperp_iw << 4*S.D0_iw['up|up']

S.G0_iw << inverse(iOmega_n + mu - V**2 * inverse(iOmega_n - eps0))

S.solve(   h_int=h_int,
           n_cycles  =20000,
           length_cycle = 10,
           n_warmup_cycles = 0,                
           measure_gt=True,
           measure_nnt=True,
           measure_chipmt=True,
           measure_nn=True
  )

if mpi.is_master_node():
  A = HDFArchive("ct_seg_single_bath_spinspin.out.h5",'w')
  A['D0_iw'] = S.D0_iw
  A['G_tau'] = S.G_tau
  A['K_tau'] = S.K_tau
  A['Jperp_tau'] = S.Jperp_tau
  A['chipm_tau'] = S.chipm_tau
  A['nn_tau'] = S.nn_tau
  print("nn =" , S.nn)
  del A
from triqs.utility.h5diff import h5diff
h5diff("ct_seg_single_bath_spinspin.out.h5", "ct_seg_single_bath_spinspin.ref.h5")
