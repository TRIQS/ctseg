from triqs.gf import *
from triqs.gf.descriptors import Function
from triqs.operators import *
from h5 import *
from numpy import matrix,array,zeros,sqrt
import triqs.utility.mpi as mpi
from triqs_ctseg import SolverCore as Solver

V = 0.5
U = 3.0
mu = U/2.
beta = 10.
e_h=0.5
alpha={}
alpha["up"]=0.25
alpha["down"]=0.5
split=.5


h_int = U*n('up',0)*n('up',1)+U*n('down',0)*n('down',1)


S = Solver(beta = beta,
           n_tau = 10000,
           gf_struct = [("up",2), ("down",2)],
           )

# initialize the Green function
for spin in ['up','down']:
 S.G0_iw[spin][0,0] = iOmega_n + mu + split/2 - V**2 * (inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h))
 S.G0_iw[spin][0,1] = -alpha[spin] * V**2 * (inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h))
 S.G0_iw[spin][1,1] = iOmega_n + mu - split/2 - V**2 * (inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h))
 S.G0_iw[spin][1,0] = -alpha[spin] * V**2 * (inverse(iOmega_n -e_h)  + inverse(iOmega_n + e_h))
S.G0_iw.invert()

S.solve(h_int=h_int,
        n_cycles  = 20000,
        length_cycle = 1,
        n_warmup_cycles = 1000,     
        measure_gt=True,
        measure_gw=True,
        measure_fw=True,
        )

if mpi.is_master_node():
 A = HDFArchive("offdiag_ctseg_offdiag_4_colors.out.h5",'w')
 A['G_iw'] = S.G_iw
 A['Sigma'] = S.Sigma_iw
 A['G_l'] = S.G_l
 del A
