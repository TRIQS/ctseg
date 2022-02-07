# PAM for two level hybridizing with a semicilcular band
#
from triqs.gf import *
from h5 import *
from triqs.gf.descriptors import *
from triqs.operators import *
from triqs.utility import mpi
from cthyb_matrix import SolverCore as Solver
from triqs.dos.dos import *
from triqs.dos.hilbert_transform import *
#from utils_LP import *
import time
import numpy

Loops =    1
Density_Required    =  2.0
Filename='PAM'

#Model parameters
Beta                =   40.0
U                   =   6.0
J                   =   0.0
eps_f               =  -2.5
CF_split            =   0.007
V0                  =   0.05
eps_V               =   2.0
t                   =   0.5

DOSFile             =   'dos.dat'
Norb=2
NFreqMatsubara=1025

#Initialization
Chemical_Potential  =   0.0
Mix                 =   1.0
load_previous       =   False

#CTQMC parameters
QMCcycles = 2000000
Length_cycle = 200
Warming_iterations =10000
fit_start =  Beta*2 
fit_interval = Beta*2 
n_tau =  10000

import triqs.utility.dichotomy as dichotomy
U                   =   U-2*J
mpi.report("U=%s"%(U))
mpi.report("J=%s"%(J))
mpi.report("Beta=%s"%(Beta))

e0=eps_f
e1=eps_f+CF_split

V1=V0*eps_V

V0sq=V0*V0
V1sq=V1*V1

mpi.report("Position of level-0: %s"%(e0))
mpi.report("Position of level-1: %s"%(e1))
mpi.report("Hybridization of level-0: %s"%(V0))
mpi.report("Hybridization of level-1: %s"%(V1))

#Omega_Max=NFreqMatsubara*numpy.pi*2.0/Beta
#mpi.report("Omega_Max: %s"%(Omega_Max))

HDFfilename = Filename+'.h5'

# Prepare_Run_Directory(DirectoryName = "SingleSite_density_%.2f"(Density_Required), OverWrite=False)

DOS1band = dos_from_file(DOSFile,single_orbital =1)
for i in range(len(DOS1band.eps)):
    mpi.report(" %10.5f  %10.5f"%(DOS1band.eps[i],DOS1band.rho[i]))
HT = HilbertTransform(DOS1band)


spin_arr=['up','down']
map =  {}
gf_struct=[ ('up',2) ,('down',2)]


#S = Solver_MultiBand_Segment(Beta = Beta, Norb = Norb, U_interact = U, GFStruct= GFStruct, map = map, J_Hund = J, Omega_Max=Omega_Max )
S = Solver(beta=Beta,gf_struct = gf_struct)

H_local = U * (N('up',0) * N('up',1)+N('down',0) * N('down',1)+N('up',0) *N('down',0)+N('up',0) *N('down',1)+N('up',1)*N('down',0)+N('up',1)*N('down',1))

Sigma=S.G0.copy()
Delta_iw=S.G0.copy()
G=S.G0.copy()

# Init Green function
#check if there are previous runs:
previous_runs = 0
previous_present = False
mu=Chemical_Potential
if mpi.is_master_node() and load_previous:
    ar = HDFArchive(HDFfilename,'a')
    if 'iterations' in ar:
        previous_present = True
        previous_runs = ar['iterations']
        Sigma <<= ar['SigmaF']
    if 'mu' in ar:
        mu = ar['mu']
    del ar

previous_runs    = mpi.bcast(previous_runs)
previous_present = mpi.bcast(previous_present)
mu = mpi.bcast(mu)
if previous_present: 
    mpi.report("Using stored data for initialisation")
    Sigma = mpi.bcast(Sigma)
else:
    #Hub-I start
    mpi.report("Hub-I start")
    for spin in ['up','down']:
        Sigma[spin][0,0] <<= Const(3*U/4)+A_Omega_Plus_B(B=U/4.0,Invert=True)*3.0*U*U/16.0
        Sigma[spin][1,1] <<= Const(3*U/4)+A_Omega_Plus_B(B=U/4.0,Invert=True)*3.0*U*U/16.0
    #save_Gf(Sigma,'Sigma_init')
    if (mpi.is_master_node()):
        ar = HDFArchive(HDFfilename,'a')
        ar['SigmaF'] = Sigma
        del ar

glist = lambda : [ GfImFreq(indices = [0], beta=Beta,n_points = NFreqMatsubara)  for i in spin_arr]
G=BlockGf(name_list = spin_arr,block_list = glist(),make_copies =False)

#glist_tau=GfImTime(indices=[0], n_points = n_tau, statistic='F', beta=Beta)
#Delta_tau=BlockGf(name_list = S.block_indices_pack_bosonic,block_list = glist_tau(),make_copies =False)
#assert 0
Gcond=G.copy()
G0=G.copy()
G1=G.copy()
G01=G.copy()

def charge(mu,G0,G1,Gcond):

    G0at=G0.copy()
    G0at_b=G0.copy()
    G1at=G0.copy()
    G1at_b=G0.copy()
    Gc0=G0.copy()

    Sig0=G0.copy()
    Sig1=G0.copy()
    Sig01=G0.copy()
    G0at_b <<=  A_Omega_Plus_B(B=mu-e0,Invert=True)
    G1at_b <<=  A_Omega_Plus_B(B=mu-e1,Invert=True)

    for spin in ['up','down']:
        Sig0[spin] <<= Sigma[spin][0,0]
        Sig1[spin] <<= Sigma[spin][1,1]
        Sig01[spin] <<= Sigma[spin][0,1]

    G0at <<= inverse( inverse(G0at_b) - Sig0)
    G1at <<= inverse( inverse(G1at_b) - Sig1)


    SigFac = G0.copy()
    SigFac <<= Const(1.0)
    SigFac <<= inverse (SigFac-Sig01*Sig01*G0at*G1at)
    SigC = (V0sq*G0at+V1sq * G1at+2.0*V0*V1*Sig01*G0at* G1at)*SigFac

    for spin in ['up','down']:
        Gcond[spin] <<= HT(SigC[spin],mu=mu)
    # recalculate S.G0
    Gc0 <<=  A_Omega_Plus_B(B=mu,Invert=False)
    for spin in ['up','down']:
        Gc0[spin] <<=  inverse(Gc0[spin] - t*t*Gcond[spin]-V0sq*G0at_b[spin]-V1sq*G1at_b[spin])
        S.G0[spin][0,0] <<= G0at_b[spin]+V0sq*G0at_b[spin]*G0at_b[spin]*Gc0[spin]
        S.G0[spin][1,1] <<= G1at_b[spin]+V1sq*G1at_b[spin]*G1at_b[spin]*Gc0[spin]
        S.G0[spin][1,0] <<= V0*V1*G0at_b[spin]*G1at_b[spin]*Gc0[spin]
        S.G0[spin][0,1] <<= S.G0[spin][1,0]

    G0 <<= V0sq+2.0*V0*V1*Sig01*G1at+V1sq*G1at*G1at*Sig01*Sig01
    G0 <<= G0at*SigFac+G0at*G0at*G0*SigFac*SigFac*Gcond
    G1 <<= V1sq+2.0*V0*V1*Sig01*G0at+V0sq*G0at*G0at*Sig01*Sig01
    G1 <<= G1at*SigFac+G1at*G1at*G1*SigFac*SigFac*Gcond

    n0=G0.total_density()
    n1=G1.total_density()
    nc=Gcond.total_density()

    charge=n0+n1+nc

    return charge


def symmSigma():
    ss= 0.5*(Sigma['up'][0,0]+Sigma['down'][0,0])
    Sigma['up'][0,0] <<= ss
    Sigma['down'][0,0] <<= ss
    ss= 0.5*(Sigma['up'][1,1]+Sigma['down'][1,1])
    Sigma['up'][1,1] <<= ss
    Sigma['down'][1,1] <<= ss
    ss= 0.25*(Sigma['up'][0,1]+Sigma['down'][0,1]+Sigma['up'][1,0]+Sigma['down'][1,0])
    Sigma['up'][0,1] <<= ss
    Sigma['up'][1,0] <<= ss
    Sigma['down'][0,1] <<= ss
    Sigma['down'][1,0] <<= ss


# DMFT loop
for IterationNumber in range(1,Loops+1) :

    itn = previous_runs + IterationNumber

    symmSigma()
  
    F = lambda mu : charge(mu,G0,G1,Gcond)

    if Density_Required :
        Chemical_potential =  dichotomy.dichotomy(function = F,
                                                x_init = mu,y_value =Density_Required,
                                                precision_on_y = 0.000001,delta_x=0.5,
                                                max_loops = 100, x_name="Chemical_Potential", y_name= "Total Density",
                                                verbosity = 3)[0]
    else:
        mpi.report("No adjustment of chemical potential\nTotal density  = %.3f"%F(Chemical_Potential))
        charge(mu,G0,G1,Gcond)

    charge(Chemical_potential,G0,G1,Gcond)
    n0=G0.total_density()
    n1=G1.total_density()
    nc=Gcond.total_density()

    mpi.report("Occupancy of level-0: %s"%(n0))
    mpi.report("Occupancy of level-1: %s"%(n1))
    mpi.report("Occupancy of conduction band: %s"%(nc))

    mpi.report("Levels Occ.: %s  %s"%(n0,n1))

#    save_Gf(S.G0,'G0imp')
    atime=time.time()
    S.solve(H_local = H_local, 
            quantum_numbers = {
              'Nup0' : N('up',0),   
              'Nup1' : N('up',1),
              'Ndown0' : N('down',0),
              'Ndown1' : N('down',1)},
            n_cycles=QMCcycles,
            length_cycle = Length_cycle,
            n_warmup_cycles =Warming_iterations,
            n_legendre = 100,                         # Number of Legendre coefficients
            random_name = 'mt19937',                 # Name of the random number generator
            use_segment_picture = False              # Use the segment picture
        )

    btime=time.time()
    mpi.report("Impurity  problem solved in %s sec"%(btime-atime))
#    save_Gf(S.G,'Gimp')
#    save_Gf(S.Delta_tau,'Delta_tau')
#    save_Gf(S.G_tau,'G_tau')

    Sigma <<= inverse(S.G0) - inverse(S.G)

    mpi.report("Orbital densities of impurity GF:")
    dm = S.G.density()
    for s,gf in S.G:
        mpi.report("Orbital %s: %s"%(s,dm[s])) 


    # Now mix Sigma with factor Mix:
    if ((IterationNumber>1) or (previous_present)):
        if (mpi.is_master_node()):
            ar = HDFArchive(HDFfilename,'a')
            mpi.report("Mixing Sigma and G with factor %s"%Mix)
            Sigma <<= Mix * Sigma + (1.0-Mix) * ar['SigmaF']
            S.G <<= Mix * S.G + (1.0-Mix) * ar['Gimp']
            del ar
        Sigma = mpi.bcast(Sigma)
        S.G = mpi.bcast(S.G)

    
    if (mpi.is_master_node()):
        ar = HDFArchive(HDFfilename,'a')
        ar['iterations'] = itn
        ar['SigmaF'] = Sigma
        ar['Gimp'] = S.G
        ar['G_tau'] = S.G_tau
        ar['G0imp'] = S.G0
        ar['Delta_tau'] = S.Delta_tau
        ar['G0']=G0
        ar['G1']=G1
        ar['Gcond']=Gcond
        ar['mu']=mu
        del ar









