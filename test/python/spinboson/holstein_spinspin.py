import numpy
from numpy import array,zeros
from numpy import sinh,cosh, cos, sin, exp, arctan
from math import sqrt, pi
from numpy import linspace

import triqs.utility.mpi as mpi
from triqs.gf import *
from triqs.gf.descriptors import Function
from triqs.operators import *
from h5 import *
from triqs_ctseg import SolverCore as Solver

def initialize(U,beta,mu,V,g,w0, S):
   ''' K_tau for a single holstein boson at energy w0 and coupling strength g'''

   sigz = lambda s1, s2 : 1 if s1==s2 else -1
   for n1 in ["up","down"]:
     for n2 in ["up","down"]:
         S.D0_iw[n1+"|"+n2] << Function(lambda w: g**2*2*w0/(w**2-w0**2))*sigz(n1,n2)
   S.Jperp_iw << 4*S.D0_iw["up|up"][0,0]

   iOm = lambda n,beta: (2*n+1)*pi/beta
   Delta_vec = 1j*V**2/2.*array([-2*arctan(1/iOm(n,beta)) for n in range(-1025,1025)])  
   Delta=S.G0_iw.copy()
   f_coeff=TailGf(1,1,0,-1)
   for n,gr in Delta:
    gr.data[:,0,0]= Delta_vec
    gr.fit_tail(f_coeff,5,10,70)

   for n,gr in S.G0_iw:  gr << inverse(iOmega_n + mu - Delta[n])

def run_spinspin_holstein(beta, V, g, time_max, U, split_regroup_v1=False, split_regroup_v2=False):

 from_scratch = True
 mu = U/2.
 w0=1.0

 S = Solver(beta = beta,                              
   gf_struct = [("up",1), ("down",1)],
            n_tau = 10000,
            n_tau_nn = 1000
            )
 Umat = numpy.matrix( [[0,U], [U,0] ])

 if from_scratch:
  initialize(U,beta,mu,V,g,w0, S)
 else: ## or load...
  A0 = HDFArchive("single_bath_seg.input.h5",'r')
  S.G0_iw = A0['G0_iw'] 
  S.D0_iw = A0['D0_iw'] 
  S.Jperp_iw = A0['Jperp_iw'] 
  del A0

 S.solve(   
            h_int = U*n('up',0)*n('down',0),
            n_cycles  =1000000000000000,
            length_cycle = 100,
            n_warmup_cycles = 10000,                
            max_time=time_max,
            move_move=False,
            move_insert_segment=True,
            move_remove_segment=True,
            move_split_spin_segment=split_regroup_v1,
            move_group_into_spin_segment=split_regroup_v1,
            move_split_spin_segment2=split_regroup_v2,
            move_group_into_spin_segment2=split_regroup_v2,
            measure_gt=True,
            measure_ft=True,
            measure_gw=False,
            measure_fw=False,
            measure_nnt=False,
            measure_chipmt=False,
            measure_nnw=False,
            measure_hist=True,
            measure_hist_composite=True,
            measure_nn=True,
   )

 if mpi.is_master_node():
   if from_scratch:
    A0 = HDFArchive("single_bath_seg.input.h5",'w')
    A0['G0_iw'] = S.G0_iw
    A0['D0_iw'] = S.D0_iw
    A0['Delta_tau'] = S.Delta_tau
    A0['K_tau'] = S.K_tau
    A0['Kperpprime_tau'] = S.Kperpprime_tau
    A0['Kprime_tau'] = S.Kprime_tau
    A0['Jperp_iw'] = S.Jperp_iw
    A0['Jperp_tau'] = S.Jperp_tau
    del A0
   A = HDFArchive("single_bath_seg.output.h5",'w')
   A['histogram'] = S.histogram
   A['histogram_composite'] = S.histogram_composite
   A['G_tau'] = S.G_tau
   A['F_tau'] = S.F_tau
   A['chipm_tau'] = S.chipm_tau
   A['nn_tau'] = S.nn_tau
   A['g0w'] = S.G0_iw
   A['sigma_w'] = S.Sigma_iw
   A['dw'] = S.D0_iw
   A['nnw'] = S.nn_iw
   A['nn'] = S.nn
   print("nn =" , S.nn)
   del A
