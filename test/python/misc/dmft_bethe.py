################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Ferrero, O. Parcollet
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

from h5 import *
from triqs.gf import *
#from triqs.dmft import DMFTLoopGeneric
import numpy
import triqs.utility.mpi as mpi
#
#  Example of DMFT single site solution with CTQMC
#
n_loops=10

# set up a few parameters
Half_Bandwidth= 1.0
U = 3.5
Chemical_Potential = U/2.0
Beta = 100.

# Construct a CTQMC solver

gf_struct = [("up",1), ("down",1)]

from ctqmc_seg import SolverCore as Solver   # imports the solver class
S = Solver(beta = Beta,
           gf_struct = gf_struct,
           n_tau = 10000,
           )

# Interaction matrix
Int = {}
Int['up'] = {}
Int['down'] = {}
Int['up']['down'] = U
Int['down']['up'] = U

Umat = numpy.zeros([2,2])
i = 0
for n1,g1 in S.G0:
 j = 0
 for n2,g2 in S.G0:
    try:
      Umat[i][j] = Int[n1][n2]
    except: pass
    j += 1
 i += 1

# init the Green function
G = S.G0.copy()
G <<= SemiCircular(Half_Bandwidth)

for IterationNumber in range(n_loops):
      print("###### LOOP ",IterationNumber)
      # Compute S.G0 with the self-consistency condition while imposing paramagnetism 
      for name, g0block in S.G0:
        g0block <<= inverse( iOmega_n + Chemical_Potential - (Half_Bandwidth/2.0)**2  * G[name] )
      
      # Run the solver
      S.solve(U=Umat,
              dynamical_U = False, 
              n_cycles  = 500000,
              length_cycle = 10,
              n_warmup_cycles = 5000,
              measure_nn = True,
              move_move = True
              )

      gt = 0.5 * ( S.G_tau['up'] + S.G_tau['down'] )
      G <<= Fourier(gt)
      
      # Some intermediate saves
      if mpi.is_master_node():
        R = HDFArchive("single_site_bethe.h5")
        R["G_tau-%s"%IterationNumber] = S.G_tau
        R["Delta_tau-%s"%IterationNumber] = S.Delta_tau
        del R
