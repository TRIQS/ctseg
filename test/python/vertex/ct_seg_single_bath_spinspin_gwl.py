from triqs.gf import *
from h5 import *
import triqs.utility.mpi as mpi
from triqs_ctseg import SolverCore as Solver
from triqs.gf.descriptors import Function
from numpy import array,zeros, sinh,cosh, matrix

filename="test_spinspin"
n_cycles=500000
n_freq=200
n_w_b_vertex = n_freq/4
n_w_f_vertex = n_freq/2
U = 1.5
beta = 20.0
mu = U/2.
eps0 = 0.0
V = 1.0
w0=0.5
l=0.0
spinboson=False

g=0.#.5  #spin-boson coupling
gperp=0.#.5
S = Solver(beta = beta, gf_struct = [("up",1), ("down",1)],
           n_tau = 10000,
           n_w_b_nn = n_freq, n_iw = n_freq)

Umat = matrix([[0.,U],[U,0.]])
S.G0 <<= inverse(iOmega_n + mu - V**2 * inverse(iOmega_n - eps0))

Kz = lambda tau : -(g/w0)**2*(cosh((beta/2-tau)*w0)-cosh(beta*w0/2))/sinh(beta*w0/2)
Kzprime = lambda tau : g**2/w0 * sinh((beta/2-tau)*w0)/sinh(beta*w0/2)
K0 = lambda tau :  (-l**2 / w0**2) * (cosh(w0*(beta/2.0-tau)) - cosh(w0*beta/2.0)) / sinh(w0*beta/2.0) 
K0prime = lambda tau : (l**2 / w0) * sinh(w0*(beta/2.0-tau)) / sinh(w0*beta/2.0)  
Jz_tau = lambda tau : -gperp**2* cosh((beta/2-tau)*w0)/sinh(beta*w0/2) 
#
#S.Kz_tau <<= Function(lambda tau: array([[1,-1],[-1,1]]) * Kz(tau) /4 )
#S.Kzprime_tau <<= Function(lambda tau: array([[1,-1],[-1,1]]) * Kzprime(tau)/4 )
#S.Jperp_tau <<= Function(lambda tau:  Jz_tau(tau) )
#S.K0_tau <<= Function(lambda tau: array([[1,1],[1,1]]) *K0(tau))
#S.K0prime_tau <<= Function(lambda tau: array([[1,1],[1,1]]) *K0prime(tau))

S.D_iw["up","up"]     <<= l**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))  - .25* g**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))
S.D_iw["down","down"] <<= l**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))  - .25* g**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))
S.D_iw["up","down"]   <<= l**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))  + .25* g**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))
S.D_iw["down","up"]   <<= l**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))  + .25* g**2*(inverse(iOmega_n - w0) - inverse(iOmega_n + w0))

S.solve(U=Umat,
        n_cycles  = n_cycles,
        length_cycle = 100,
        n_warmup_cycles = 500,
        n_w_f_vertex = n_w_f_vertex,
        n_w_b_vertex = n_w_b_vertex,
#random_seed = 34788+928374*rank,
        evaluate_vertex = True,
        measure_gw = True,
        measure_gt = True,
        measure_nnw = True,
        measure_g2w = True,
        measure_f2w = False,
        measure_fw = True,
        measure_nn=True,
        move_move=False,
        move_swap_empty_lines=False,        
        )

if mpi.is_master_node():
  A = HDFArchive(filename+".out.h5",'w')
  A['gw'] = S.G_iw
  A['gt'] = S.G_tau
  A['g0w'] = S.G0
  A['fw'] = S.F_iw
  A['dw'] = S.D_iw
  A['nnw'] = S.nn_iw
  A['nn'] = S.nn
  A['Sw'] = S.Sigma_iw
  A['U_mat'] = Umat
#  A['Sw_from_bubble'] = S.Sigma_from_bubble_iw
  print(S.nn)

S.print_3leg_vertices(filename+".vertex.out.h5", 0)

 
from triqs.utility.h5diff import h5diff
h5diff("ct_seg_single_bath_spinspin_gwl.out.h5", "ct_seg_single_bath_spinspin_gwl.ref.h5")
