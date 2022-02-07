# PAM for two level hybridizing with a semicircular band
#
from triqs.gf import *
from h5 import *
from triqs.gf.descriptors import *
from triqs.utility import mpi
from triqs_ctseg import SolverCore as Solver
from triqs.dos.dos import *
from triqs.dos.hilbert_transform import *
import time
from numpy import *

Loops =    1
Density_Required    =  2.0
Filename='PAM'

#Model parameters
Beta                =   900.0
U                   =   6.0
J                   =   0.0
eps_f               =  -2.5
CF_split            =   0.007
V0                  =   0.03
eps_V               =   2.0
t                   =   0.5

DOSFile             =   'dos.dat'
Norb=2
NFreqMatsubara=1025

#Initialization
Chemical_potential  =   0.0
Mix                 =   1.0
load_previous       =   False

#CTQMC parameters
QMCcycles = 2350000
Length_cycle = 200
Warming_iterations =10000
measure_gt=False
measure_gw= not measure_gt
improved_estimator=not measure_gt
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

HDFfilename = Filename+'.h5'

DOS1band = dos_from_file(DOSFile,single_orbital =1)
#for i in range(len(DOS1band.eps)):
#    mpi.report(" %10.5f  %10.5f"%(DOS1band.eps[i],DOS1band.rho[i]))
HT = HilbertTransform(DOS1band)

spin_arr=['up','down']
map =  {}

gf_struct = [('up',2),('down',2)]

S = Solver(beta=Beta,n_tau=n_tau, gf_struct = gf_struct)

Umat=numpy.zeros([4,4])
Umat.fill(U)
for i in range(4):
    Umat[i,i]=0.0

if mpi.is_master_node():
    print('Umat: ')
    print(Umat)
    print(Beta)
    print("total iterations: "+str(QMCcycles*Length_cycle))

Sigma=S.G0_iw.copy()
Delta_iw=S.G0_iw.copy()
G=S.G0_iw.copy()

#set_coupling()
# Init Green function
#check if there are previous runs:
previous_runs = 0
previous_present = False
mu=Chemical_potential
if mpi.is_master_node() and load_previous:
    ar = HDFArchive(HDFfilename,'a')
    if 'iterations' in ar:
        previous_present = True
        previous_runs = ar['iterations']
        Sigma << ar['SigmaF']
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
        Sigma[spin][0,0] << Const(3*U/4)+A_Omega_Plus_B(B=U/4.0,Invert=True)*3.0*U*U/16.0
        Sigma[spin][1,1] << Const(3*U/4)+A_Omega_Plus_B(B=U/4.0,Invert=True)*3.0*U*U/16.0
    if (mpi.is_master_node()):
        ar = HDFArchive(HDFfilename,'a')
        ar['SigmaF'] = Sigma
        del ar


glist = lambda : [ GfImFreq(indices = [0], beta=Beta,n_points = NFreqMatsubara)  for i in spin_arr]
Gcond=BlockGf(name_list = spin_arr,block_list = glist(),make_copies =False)

G0=Gcond.copy()
G1=Gcond.copy()

def charge(mu,G0,G1,Gcond):

    G0at=G0.copy()
    G0at_b=G0.copy()
    G1at=G0.copy()
    G1at_b=G0.copy()
    Gc0=G0.copy()

    Sig0=G0.copy()
    Sig1=G0.copy()
    Sig01=G0.copy()
    G0at_b <<  A_Omega_Plus_B(B=mu-e0,Invert=True)  #(Omega * A + B).inverse
    G1at_b <<  A_Omega_Plus_B(B=mu-e1,Invert=True)


    for spin in ['up','down']:
        Sig0[spin] << Sigma[spin][0,0]
        Sig1[spin] << Sigma[spin][1,1]
        Sig01[spin] << Sigma[spin][0,1]


    G0at << inverse( inverse(G0at_b) - Sig0)  #S_gamma in the paper
    G1at << inverse( inverse(G1at_b) - Sig1)

    SigFac = G0.copy()
    SigFac << Const(1.0)
    SigFac << inverse (SigFac-Sig01*Sig01*G0at*G1at)
    SigC = (V0sq*G0at+V1sq * G1at+2.0*V0*V1*Sig01*G0at* G1at)*SigFac

    for spin in ['up','down']:
        Gcond[spin] << HT(SigC[spin],mu=mu)

    # recalculate S.G0_iw
    Gc0 <<  A_Omega_Plus_B(B=mu,Invert=False)
    for spin in ['up','down']:
        Gc0[spin] <<  inverse(Gc0[spin] - t*t*Gcond[spin]-V0sq*G0at_b[spin]-V1sq*G1at_b[spin])
        S.G0_iw[spin][0,0] << G0at_b[spin]+V0sq*G0at_b[spin]*G0at_b[spin]*Gc0[spin]
        S.G0_iw[spin][1,1] << G1at_b[spin]+V1sq*G1at_b[spin]*G1at_b[spin]*Gc0[spin]
        S.G0_iw[spin][1,0] << V0*V1*G0at_b[spin]*G1at_b[spin]*Gc0[spin]
        S.G0_iw[spin][0,1] << S.G0_iw[spin][1,0]


    G0 << V0sq+2.0*V0*V1*Sig01*G1at+V1sq*G1at*G1at*Sig01*Sig01
    G0 << G0at*SigFac+G0at*G0at*G0*SigFac*SigFac*Gcond
    G1 << V1sq+2.0*V0*V1*Sig01*G0at+V0sq*G0at*G0at*Sig01*Sig01
    G1 << G1at*SigFac+G1at*G1at*G1*SigFac*SigFac*Gcond


    n0=G0.total_density()
    n1=G1.total_density()
    nc=Gcond.total_density()

    charge=n0+n1+nc

    return charge


def symmSigma():
    ss= 0.5*(Sigma['up'][0,0]+Sigma['down'][0,0])
    Sigma['up'][0,0] << ss
    Sigma['down'][0,0] << ss
    ss= 0.5*(Sigma['up'][1,1]+Sigma['down'][1,1])
    Sigma['up'][1,1] << ss
    Sigma['down'][1,1] << ss
    ss= 0.25*(Sigma['up'][0,1]+Sigma['down'][0,1]+Sigma['up'][1,0]+Sigma['down'][1,0])
    Sigma['up'][0,1] << ss
    Sigma['up'][1,0] << ss
    Sigma['down'][1,0] << ss
    Sigma['down'][0,1] << ss

# DMFT loop
for IterationNumber in range(1,Loops+1) :

    itn = previous_runs + IterationNumber

    symmSigma()

    F = lambda mu : charge(mu,G0,G1,Gcond)

    if Density_Required :
        Chemical_potential =  dichotomy.dichotomy(function = F,
                                                x_init = mu,y_value =Density_Required,
                                                precision_on_y = 0.000001,delta_x=0.5,
                                                max_loops = 100, x_name="Chemical_potential", y_name= "Total Density",
                                                verbosity = 3)[0]
    else:
        mpi.report("No adjustment of chemical potential\nTotal density  = %.3f"%F(Chemical_potential))
        charge(mu,G0,G1,Gcond)

    charge(Chemical_potential,G0,G1,Gcond)
    n0=G0.total_density()
    n1=G1.total_density()
    nc=Gcond.total_density()

    mpi.report("\nIteration %s\n"%(itn))
    mpi.report("Occupancy of level-0: %s"%(n0))
    mpi.report("Occupancy of level-1: %s"%(n1))
    mpi.report("Occupancy of conduction band: %s"%(nc))

    mpi.report("Levels Occ.: %s  %s"%(n0,n1))


    atime=time.time()
    try:
     S.solve(U=Umat, n_cycles=QMCcycles,length_cycle = Length_cycle,
            n_warmup_cycles =Warming_iterations,
            measure_gw=False,measure_gt=measure_gt,measure_fw=False,
            #measure_gw=measure_gw,measure_gt=measure_gt,measure_fw=improved_estimator,
            move_move=True,move_swap_empty_lines=True
            )
    except:
       print("error")
       if (mpi.is_master_node()):
        ar = HDFArchive(HDFfilename,'a')
        ar['Delta_tau'] = S.Delta_tau
        del ar

    btime=time.time()
    mpi.report("Impurity  problem solved in %s sec"%(btime-atime))
    if (mpi.is_master_node()):
        ar = HDFArchive(HDFfilename,'a')
        ar['Gimp'] = S.G_iw
        ar['G_tau'] = S.G_tau
        del ar


    # get G,Sigma in omega
    if not improved_estimator :
        for s,s1,g in S.G_tau:
            g.tail.zero()
            if s==s1:
                g.tail[1] = numpy.identity(g.N1)
            G[s,s1].set_from_fourier(g)
            if measure_gt : G[s,s1].set_from_fourier(g)

        Sigma << inverse(S.G0_iw) - inverse(G)
        save_Gf(Sigma,'Sigma')
        fixed_coeff = TailGf(1,1,1,0);
        fixed_coeff[0] = array([[0.]]);
        fit_n_moments=4;
        for n,g in Sigma:
           g.fit_tail(fixed_coeff, fit_n_moments, fit_start, fit_start+fit_interval);

    else:
        
        Sigma << S.Sigma_iw
        fixed_coeff = TailGf(2,2,1,0);
        fixed_coeff[0] = array([[0.,0.],[0.,0.]]);
        fit_n_moments=4;
        #for n,g in Sigma:
        #    g.fit_tail(fixed_coeff, fit_n_moments, int(fit_start), int(fit_start+fit_interval)); 
 
    G=S.G_iw

    mpi.report("Orbital densities of impurity GF:")
    dm = G.density()
    if (mpi.is_master_node()): print(dm)
    #for s,s1,gf in G:
    #    if s==s1 :
    #        mpi.report("Orbital %s: %s"%(s,dm[s])) 


    # Now mix Sigma with factor Mix:
    if (mpi.is_master_node()):
        ar = HDFArchive(HDFfilename,'a')
        mpi.report("Mixing Sigma and G with factor %s"%Mix)
        Sigma << Mix * Sigma + (1.0-Mix) * ar['SigmaF']
        del ar
    Sigma = mpi.bcast(Sigma)

    
    if (mpi.is_master_node()):
        ar = HDFArchive(HDFfilename,'a')
        ar['Delta_tau'] = S.Delta_tau
        ar['G0_iw'] = S.G0_iw
        ar['iterations'] = itn
#        ar['SigmaF'] = Sigma
        ar['Sigma_iw'] = S.Sigma_iw
        ar['Gimp'] = G
        ar['G0']=G0
        ar['G1']=G1
        ar['Gcond']=Gcond
        ar['mu']=mu
        del ar




