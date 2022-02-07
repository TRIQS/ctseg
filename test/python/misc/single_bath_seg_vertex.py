from triqs.gf import *
from triqs.gf.descriptors import Function
from h5 import *
import numpy
import triqs.utility.mpi as mpi

# set up a few parameters
U = 2.0
Beta = 10.0
mu = 1.0
eps0 = 0.0

# Construct a Segment solver
from triqs_ctseg import SolverCore as Solver
S = Solver(beta = Beta,
           gf_struct = [("up",1), ("down",1)],
           n_tau= 10000,
           n_legendre_g = 20,
           n_tau_nn = 1000,
           n_w_b = 2,
           n_w_f_vertex = 2
           )

# Interaction matrix
Umat = numpy.zeros([2,2])
for (i,(n1,g1)) in enumerate (S.G0) : 
  for (j,(n2,g2)) in enumerate (S.G0) :
      Umat[i][j] = U if set((n1,n2)) == set(['down','up']) else 0

l=1.0
omega0=1.0

#Setting S.K_tau and S.Kprime
S.K_tau <<= Function(lambda tau: numpy.array([[1,1],[1,1]]) * (-l**2 / omega0**2) * (numpy.cosh(omega0*(Beta/2.0-tau)) - numpy.cosh(omega0*Beta/2.0)) / numpy.sinh(omega0*Beta/2.0) )
S.Kchprime_tau <<= Function(lambda tau: numpy.array([[1,1],[1,1]]) * (l**2 / omega0) * numpy.sinh(omega0*(Beta/2.0-tau)) / numpy.sinh(omega0*Beta/2.0) )

# init the Green function
S.G0 <<= inverse(iOmega_n + mu - inverse(iOmega_n - eps0))

numpy.savetxt('Ktau_.dat', S.K_tau.data[:,0,0])

S.solve(U=Umat,
        dynamical_U=False,
        improved_estimator=0,
        n_cycles  = 1000000,
        length_cycle = 10,
        max_time = 10,
        n_warmup_cycles = 1000,
        measure_nn=True,
        measure_nnt=False,
        measure_nnw=True,
        measure_gw=False,
        measure_g2w=True,
        measure_g3w=True,
        measure_gl=False)

if mpi.is_master_node():
  A = HDFArchive("single_bath_seg.output.h5",'w')
  A['G_tau'] = S.G_tau

  numpy.savetxt('gtau_up.dat', S.G_tau['up'].data)
  numpy.savetxt('gtau_dn.dat', S.G_tau['down'].data)

  numpy.savetxt('delta_up.dat', S.Delta_tau['up'].data)
  numpy.savetxt('delta_dn.dat', S.Delta_tau['down'].data)

  numpy.savetxt('ftau_up.dat', S.F_tau['up'].data)
  numpy.savetxt('ftau_dn.dat', S.F_tau['down'].data)

  numpy.savetxt('gl_up.dat', S.G_legendre['up'].data)
  numpy.savetxt('gl_dn.dat', S.G_legendre['down'].data)

  numpy.savetxt('fl_up.dat', S.F_legendre['up'].data)
  numpy.savetxt('fl_dn.dat', S.F_legendre['down'].data)

  numpy.savetxt('gw_up.dat', S.G_iw['up'].data.imag)
  numpy.savetxt('gw_dn.dat', S.G_iw['down'].data.imag)
  numpy.savetxt('gw_up_r.dat', S.G_iw['up'].data.real)
  numpy.savetxt('gw_dn_r.dat', S.G_iw['down'].data.real)

  numpy.savetxt('fw_up.dat', S.F_iw['up'].data.imag)
  numpy.savetxt('fw_dn.dat', S.F_iw['down'].data.imag)
  numpy.savetxt('fw_up_r.dat', S.F_iw['up'].data.real)
  numpy.savetxt('fw_dn_r.dat', S.F_iw['down'].data.real)

  numpy.savetxt('sw_up.dat', S.Sigma_iw['up'].data.imag)
  numpy.savetxt('sw_dn.dat', S.Sigma_iw['down'].data.imag)
  numpy.savetxt('sw_up_r.dat', S.Sigma_iw['up'].data.real)
  numpy.savetxt('sw_dn_r.dat', S.Sigma_iw['down'].data.real)

  numpy.savetxt('nntau_00.dat', S.nn_tau.data[:,0,0])
  numpy.savetxt('nntau_01.dat', S.nn_tau.data[:,0,1])
  numpy.savetxt('nntau_10.dat', S.nn_tau.data[:,1,0])
  numpy.savetxt('nntau_11.dat', S.nn_tau.data[:,1,1])

  numpy.savetxt('nnw_00.dat', numpy.real(S.nn_iw.data[:,0,0]))
  numpy.savetxt('nnw_01.dat', numpy.real(S.nn_iw.data[:,0,1]))
  numpy.savetxt('nnw_10.dat', numpy.real(S.nn_iw.data[:,1,0]))
  numpy.savetxt('nnw_11.dat', numpy.real(S.nn_iw.data[:,1,1]))

  numpy.savetxt('nnw_00_i.dat', numpy.imag(S.nn_iw.data[:,0,0]))
  numpy.savetxt('nnw_01_i.dat', numpy.imag(S.nn_iw.data[:,0,1]))
  numpy.savetxt('nnw_10_i.dat', numpy.imag(S.nn_iw.data[:,1,0]))
  numpy.savetxt('nnw_11_i.dat', numpy.imag(S.nn_iw.data[:,1,1]))


  numpy.savetxt('kt.dat', S.K_tau.data[:,0,0])
  numpy.savetxt('kprimet.dat', S.Kchprime_tau.data[:,0,0])

  print(S.nn)
