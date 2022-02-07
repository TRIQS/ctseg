from triqs.gf import *
from triqs.operators import *
from triqs.gf.descriptors import Function
from h5 import *
from numpy import matrix,array,zeros,sqrt
from triqs_ctseg import SolverCore as Solver
import triqs.utility.mpi as mpi

V = 0.5
U = 3.0
mu = U/2.
beta = 10.
e_h=0.5
alpha=0.5
split=.05


S = Solver(beta = beta,
           n_tau = 10000,
           gf_struct = [('up',2)],
           )

h_int = U*n('up',0)*n('up',1)

# initialize the Green function
S.G0_iw['up'][0,0] = iOmega_n + mu + split/2 - V**2 * (inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h)) 
S.G0_iw['up'][0,1] =  - alpha * V**2 * (inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h))
S.G0_iw['up'][1,1] = iOmega_n + mu - split/2 - V**2 * (inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h)) 
S.G0_iw['up'][1,0] =  - alpha * V**2 * (inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h))
S.G0_iw.invert()

S.solve(h_int=h_int,
        n_cycles  = 5000,
        length_cycle = 100,
        n_warmup_cycles = 1000,
        measure_gt=True,
        measure_gw=True,
        measure_fw=True,
        measure_sign=True,
        )

if mpi.is_master_node():
 A = HDFArchive("offdiag_ctseg.out.h5",'w')
 A['G_l'] = S.G_l
 A['G_iw'] = S.G_iw
 A['Sigma'] = S.Sigma_iw
 A['sign'] = S.average_sign
 del A

from triqs.utility.h5diff import h5diff
h5diff("offdiag_ctseg.out.h5", "offdiag_ctseg.ref.h5")
