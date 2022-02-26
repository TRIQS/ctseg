from triqs.gf import *
from triqs.gf.descriptors import Function
from triqs.operators import *
from h5 import *
from triqs_ctseg import SolverCore
import triqs.utility.mpi as mpi
from numpy import array,zeros,sinh, exp, arctan, linspace
from math import sqrt, pi

def initialize(S, n_w, beta, mu, V, g, wc, s, l, w0, sign=1.0):
   ''' 
       S impurity solver instance
       n_w nb of matsubara frequencies
       beta inverse temperature
       mu chemical potential
       V hybridization strength

       g coupling to the spin subohmic bath
       wc cutoff frequency for subohmic bath
       s exponent of the subohmic bath

       l charge
       w0 holstein frequency for charge-charge interactions

       sign sign in front of expression for spin-spin interactions

       Note:
       Kz defined as:
       Kz''(tau) = Dz(tau) ; Kz(0^{+}) = Kz(beta^{-}) = 0
       J_\perp(\tau) = 4 Dz(\tau) in the rotation invariant case
   '''

   nu = lambda n,beta : (2*n*pi/beta)
   Eps=linspace(0,wc,5000)
   Dz_expr = lambda w : sign*(-1.)*(.5*g**2) *((s+1)*wc**(-s-1))*sum([eps**(s+1)/(w.imag*w.imag+eps**2) for eps in Eps])*(Eps[1]-Eps[0])
   Dz0 = sign*(-1.)*(.5*g**2)*(s+1.)/s/wc #Dz_expr ill-defined at w =0 
   Dz_expr_reg = lambda w : Dz_expr(w) if abs(w)>0.0 else Dz0
   Dz_iw = GfImFreq(indices=['perp'], beta=beta,n_points = n_w, statistic='Boson')
   Dz_iw << Function(lambda w : Dz_expr_reg(w))

   S.Jperp_iw <<  4*Dz_iw

   D0_expr = lambda w : l**2*(1./(w - w0) - 1./(w + w0)) 

   S.D0_iw['up|up'] << Function(lambda w : D0_expr(w) + Dz_expr_reg(w))
   S.D0_iw['down|down'] << Function(lambda w : D0_expr(w) + Dz_expr_reg(w))
   S.D0_iw['up|down'] << Function(lambda w : D0_expr(w) - Dz_expr_reg(w))
   S.D0_iw['down|up'] << Function(lambda w : D0_expr(w) - Dz_expr_reg(w))

   iOm = lambda n,beta: (2*n+1)*pi/beta
   Delta_vec = 1j*V**2/2.*array([-2*arctan(1/iOm(n,beta)) for n in range(-n_w,n_w)])  
   Delta=S.G0_iw.copy()
   for n,gr in Delta:
    gr.data[:,0,0]= Delta_vec

   K0prime_expr = lambda tau : 0.0 if l==0.0 else (l**2/w0) * sinh(w0*(beta/2.0-tau)) / sinh(w0*beta/2.0) 
   for n,gr in S.G0_iw:
     gr << inverse(iOmega_n + mu - 2*K0prime_expr(0.0)  - Delta[n])


def run_spinspin(beta, Vsq, g, time_max, rank=-1, from_scratch=True, U=0.4, l=0.0, measure_fw=False, full_spin_rot_inv=True, split_regroup_v1=False,split_regroup_v2=False, sign=1.0, mu=0.2, move_move=False, move_swap_empty_lines=True, measure_gw=False, measure_g2w=False, n_w=400, filename="single_bath_seg", n_cycles=1000000000000):
 ''' 
     This function launches a run for an Anderson model with dynamical charge-charge (parametrized by l) and spin-spin interactions (subohmic bath parametrized by g),
     U is the strength of the static interactions
     Vsq the strength of the hybridization (squared)
     beta the inverse temperature  

     Note: in the rotation invariant case, one should have:
      <s_+(tau) s_-(0)>/2 = <s_z(tau) s_z(0)>
      with s_z = 0.5*(n_up - n_down)
 '''

 if Vsq==0.0:
    move_insert_segment=False
    move_remove_segment=False
    Vsq=.001
 else:
    move_insert_segment=True
    move_remove_segment=True

 V = sqrt(Vsq)
 wc=1.0
 s=0.2
 w0=0.5
 S = SolverCore(beta = beta,
		   gf_struct = [["up",1], ["down",1]],
		   n_tau= 10000,
     n_tau_nn = 10000,
		   n_w_b_nn = n_w,
		   n_iw = n_w
		   )

 if from_scratch:
  initialize(S, n_w, beta, mu, V, g, wc, s, l, w0, sign)
  if mpi.is_master_node():
   A0 = HDFArchive(filename+".input.h5", 'a') 
   A0['G0_iw'] = S.G0_iw
   A0['D0_iw'] = S.D0_iw
   A0['Jperp_iw'] = S.Jperp_iw
   del A0
 else: ## or load...
  A0 = HDFArchive(filename+".input.h5",'r')
  S.G0_iw << A0['G0_iw'] 
  S.D0_iw <<  A0['D0_iw'] 
  S.Jperp_iw << A0['Jperp_iw'] 
  del A0
 try:
  S.solve(h_int = U*n('up',0,)*n('down',0),
     n_w_f_vertex = 50,
 	   n_w_b_vertex = 50,
		   n_cycles  = n_cycles,
		   length_cycle = 100,
		   n_warmup_cycles = 10000,                
		   max_time=time_max,
     keep_Jperp_negative=False,
		   move_move=move_move,
     move_swap_empty_lines=move_swap_empty_lines,
     move_insert_segment = move_insert_segment,
		   move_remove_segment = move_remove_segment,
     move_split_spin_segment=split_regroup_v1,
     move_group_into_spin_segment=split_regroup_v1,
     move_split_spin_segment2=split_regroup_v2,
     move_group_into_spin_segment2=split_regroup_v2,
		   measure_gt=True,
		   measure_ft=True,
		   measure_gw=measure_gw,
		   measure_fw=measure_fw,
		   measure_g2w=measure_g2w,
		   measure_nnt=True,
		   measure_chipmt=True,
		   measure_nnw=True,
		   measure_hist=True,
		   measure_hist_composite=True,
		   measure_nn=True,
	  )
 except Exception as e:
   print("Error:",e)
   if from_scratch and mpi.is_master_node():
    A0 = HDFArchive(filename+".input.h5",'a')
    A0['K_tau'] = S.K_tau
    A0['Kprime_tau'] = S.Kprime_tau
    A0['Kperpprime_tau'] = S.Kperpprime_tau
    A0['Jperp_tau'] = S.Jperp_tau
    A0['Delta_tau'] = S.Delta_tau
    del A0
   raise

 if mpi.is_master_node():
   if from_scratch:
    A0 = HDFArchive(filename+".input.h5",'a')
    A0['Delta_tau'] = S.Delta_tau
    A0['K_tau'] = S.K_tau
    A0['Kprime_tau'] = S.Kprime_tau
    A0['Kperpprime_tau'] = S.Kperpprime_tau
    A0['Jperp_tau'] = S.Jperp_tau
    del A0
   A = HDFArchive(filename+".out.h5",'a')
   A['histogram'] = S.histogram
   A['histogram_composite'] = S.histogram_composite
   A['G_tau'] = S.G_tau
   A['F_tau'] = S.F_tau
   A['gw'] = S.G_iw
   A['fw'] = S.F_iw
   A['g0w'] = S.G0_iw
   A['sigma_w'] = S.Sigma_iw
   A['dw'] = S.D0_iw
   A['chipm_tau'] = S.chipm_tau
   A['nn_tau'] = S.nn_tau
   A['nnw'] = S.nn_iw
   A['nn'] = S.nn
   A['U']=U
   del A
   print("nn =" , S.nn)
   if measure_g2w:
    A = HDFArchive(filename+".g2w.h5",'a')
    A["g2w"] = S.G_2w
    del A
