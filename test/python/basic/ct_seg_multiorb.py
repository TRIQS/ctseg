#########################################################################
#### a multi orbital test for the segment-picture hyb-expansion solver ##
#########################################################################
from triqs.gf import *
from triqs.gf.descriptors import Function
from triqs.utility import mpi
from triqs.operators import *
from h5 import *
from numpy import matrix,array,zeros,sqrt

#interaction
U    = 2.0
Up   = 1.0
J    = 0.0
#Hamiltonian for the 'U-3J' model
h_int = U*(n('up1',0)*n('down1',0)+n('up2',0)*n('down2',0))+\
        (Up-J)*(n('up1',0)*n('up2',0)+n('down1',0)*n('down2',0))+\
         Up*(n('up1',0)*n('down2',0)+n('down1',0)*n('up2',0))

beta = 2.0 #inverse temperature
mu   = 0.2 #chemical potential

#hybridization levels
eps ={} 
eps["up1"] = 0.3
eps["down1"] = 0.3
eps["up2"] = 0.3
eps["down2"] = 0.3

#hybridization strengths
V ={}
V["up1"] = 1.0
V["down1"] = 1.0
V["up2"] = 1./sqrt(2.0)
V["down2"] = 1./sqrt(2.0)

# Construct a Segment solver
from triqs_ctseg import SolverCore as Solver

S = Solver(beta = beta,
           n_tau = 10000,
           gf_struct= [("up1",1), ("down1",1), ("up2",1), ("down2",1)]
           )

# initialize the Green function
for name, block in S.G0_iw:
 block << inverse(iOmega_n + mu - V[name]**2 * inverse(iOmega_n - eps[name]))

S.solve(h_int=h_int, 
        n_cycles  = 5000000,
        length_cycle = 1,
        n_warmup_cycles = 0,     
        measure_gt=True,
        measure_nn=True,
        )

A = HDFArchive("ct_seg_multiorb.out.h5",'w')
A['G_tau'] = S.G_tau
from triqs.utility.h5diff import h5diff
h5diff("ct_seg_multiorb.out.h5", "ct_seg_multiorb.ref.h5")
