from triqs.gf import *
from h5 import *
from numpy import matrix
import triqs.utility.mpi as mpi
#from triqs_ctseg import SolverCore as Solver
from triqs_ctseg import SolverCore as Solver
from triqs.gf.descriptors import Function
from numpy import array,zeros
from numpy import sinh,cosh

U = 1.0
Beta = 20.0
mu = .5
eps0 = 0.0
V = 1.0
w0=1.0
l=0.45

S = Solver(beta = Beta,   
           gf_struct = [("up",1), ("down",1)],
           n_tau = 10000,
        n_w_b_nn = 50,
        n_iw = 200,
           )

Umat = matrix([[0.,U],[U,0.]])

S.G0 <<= inverse(iOmega_n + mu - V**2 * inverse(iOmega_n - eps0))
S.K_tau <<= Function(lambda tau: array([[1,1],[1,1]]) * (-l**2 / w0**2) * (cosh(w0*(Beta/2.0-tau)) - cosh(w0*Beta/2.0)) / sinh(w0*Beta/2.0) )
S.Kchprime_tau <<= Function(lambda tau: array([[1,1],[1,1]]) * (l**2 / w0) * sinh(w0*(Beta/2.0-tau)) / sinh(w0*Beta/2.0) )
#S.D_iw['up','up'] <<= l**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))
#S.D_iw['down','up'] <<= l**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))
#S.D_iw['down','down'] <<= l**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))
#S.D_iw['up','down'] <<= l**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))

S.solve(U=Umat,
        move_move=False,
        move_swap_empty_lines=False,        
        n_cycles  = 100000,
#        n_w_f_vertex = 100,
#        n_w_b_vertex = 50,
        n_w_f_vertex = 20,
        n_w_b_vertex = 10,
        length_cycle = 100,
        n_warmup_cycles = 100,
        measure_gw = True,
        measure_gt = False,
        measure_nnw = True,
#measure_g2w = True,
        measure_g3w = True,
        evaluate_vertex = True,
        measure_nn=True)

if mpi.is_master_node():
  A = HDFArchive("ct_seg_single_bath.out.h5",'w')
  A['G'] = S.G_iw
  A['nn_w'] = S.nn_iw
  A['Sw'] = S.Sigma_iw
#A['Sw_from_bubble'] = S.Sigma_from_bubble_iw
 
from triqs.utility.h5diff import h5diff
h5diff("ct_seg_single_bath_gwl.out.h5", "ct_seg_single_bath_gwl.ref.h5")
