from triqs.dft.sumk_lda import *
from triqs.dft.converters import *
from triqs.dft.solver_multiband_seg import *
from triqs.base.gf_local import *
import time, numpy


#=====================================================
#Basic input parameters:
LDAFilename='SrMnO3'
U = 4.00
J = 0.80
Beta =10
Chemical_Pot_init = 0.0
Loops =  1                       # Number of DMFT sc-loops
Mix = 0.7                        # Mixing factor of Sigma after solution of the AIM
DC_type = 0                      # DC type: 0 FLL, 1 Held, 2 AMF
useBlocs = True                  # use bloc structure from LDA input
useMatrix = True                 # True: Slater parameters, False: Kanamori parameters U+2J, U, U-J
use_spinflip = False             # use the full rotational invariant interaction?
hfield = 0.000
prec_mu = 0.000001

QMCcycles = 25000
Length_cycle = 50
Warming_iterations = 1000
fit_start = 800
fit_interval = 200

#=====================================================

# Convert DMFT input:
# Can be commented after the first run
Converter = Wien2kConverter(filename=LDAFilename,repacking=True)
Converter.convert_dmft_input()
mpi.barrier()

#check if there are previous runs:
previous_runs = 0
previous_present = False

# if previous runs are present, no need for recalculating the bloc structure
# It has to be commented, if you run this script for the first time, starting
# from a converted h5 archive.

calc_blocs = useBlocs and (not previous_present)

# Init the SumK class
SK=SumkLDA(hdf_file=LDAFilename+'.h5',use_lda_blocks=calc_blocs,h_field = hfield)


Norb = SK.corr_shells[0][3]
l = SK.corr_shells[0][2]

# Init the Solver:
S=SolverMultiBandSeg(Beta=Beta,U_interact=U,J_Hund=J,Norb=Norb,useMatrix=useMatrix, T=SK.T[0] ,GFStruct=SK.gf_struct_solver[0],
                   map=SK.map[0], l=l, deg_orbs=SK.deg_shells[0], Omega_Max = 1025 * 2 * numpy.pi / Beta) # ,dimreps=dimreps,irep=irep)

S.N_Cycles  = QMCcycles
S.Length_Cycle = Length_cycle
S.N_Warmup_Cycles = Warming_iterations
#S.Time_Accumulation = True  
#S.Legendre_Accumulation = False
#S.N_Legendre_Coeffs = 35
S.Fit_Tails = True
S.Fitting_Frequency_Start = fit_start
S.N_Frequencies_Accumulated = fit_start + fit_interval


# DMFT loop:
for IterationNumber in range(1,Loops+1) :

      SK.symm_deg_gf(S.Sigma,orb=0)                           # symmetrise Sigma
      
      SK.put_Sigma(Sigma_imp = [ S.Sigma ])                    # put Sigma into the SumK class:

      SK.chemical_potential = -0.038061

      S.G <<= SK.extract_G_loc()[0]                            # calculation of the local Green function

      dm = S.G.density()
      mpi.report("Orbital densities of impurity Green function before solver:")
      for s,gf in S.G:
          mpi.report("Orbital %s: %s"%(s,dm[s]))

      if ( ( not previous_present ) and IterationNumber == 1):
          # Set initial double counting:
          dm = S.G.density()
          SK.set_dc( dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = DC_type)
          S.Sigma <<= SK.dc_imp[0]['up'][0,0]

      mpi.report("Total charge of Gloc : %.6f"%S.G.total_density())

      # now calculate new G0:
      S.G0 <<= inverse ( S.Sigma + inverse(S.G) )
      # Solve the impurity problem:
      atime=time.time()
      S.Solve()
      btime=time.time()
      mpi.report("Impurity  problem solved in %s sec"%(btime-atime))

      # solution done, do the post-processing:
      mpi.report("Total charge of impurity problem : %.6f"%S.G.total_density())

      # Write the final Sigma and G to the hdf5 archive:
      if (mpi.is_master_node()):
          ar = HDFArchive(LDAFilename+'.output.h5','a')
          ar['G_tau'] = S.G_tau
          del ar
