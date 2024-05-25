# Two impurity sites coupled to two bath sites.
# Delta(tau) has off-diagonal components. Can be compared to ED reference dimer_pyed.ref.h5
# Comparison to CTSEG reference is disabled for this test because results can have large 
# variations across platforms for short runs. 

from triqs.gf import *
from triqs.gf.tools import *
from triqs.gf.gf_factories import make_gf_from_fourier
from triqs.operators.util import U_matrix_kanamori, h_int_density
import h5
import triqs.utility.mpi as mpi
from triqs.utility.h5diff import h5diff

from triqs_ctseg import Solver

import numpy as np
from numpy import linalg 

# Numerical values
beta = 10.                      # Inverse temperature
mu = 0.0                        # Chemical potential
eps = np.array([0.0, 0.1])         # Impurity site energies
eps_bath = np.array([0.27, -0.4])  # Bath site energies
U = 1.                          # Density-density interaction
J = 0.2                         # Hund's coupling
V = 1                           # Bath-impurity coupling
block_names = ['up', 'dn']
n_orb = len(eps)
n_orb_bath = len(eps_bath)
n_tau = 10001
gf_struct = [(s, n_orb) for s in block_names] # Green's function structure

########################################################################
# --------------- Construct Delta(tau) and chemical_potential ---------------

# Non-interacting impurity Hamiltonian in matrix representation
h_0_mat = np.diag(eps - mu)

# Bath Hamiltonian in matrix representation
h_bath_mat = np.diag(eps_bath)

# Coupling matrix
V_mat = np.matrix([[V, V],
                [V, V]])

# Total non-interacting Hamiltonian 
h_tot_mat = np.block([[h_0_mat, V_mat],
                   [V_mat.H, h_bath_mat]])

# Non-interacting impurity Green's function
n_iw = int(10 * beta)
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
for bl, iw in product(block_names, iw_mesh):
    G0_iw[bl][iw] = linalg.inv(iw.value * np.eye(2*n_orb) - h_tot_mat)[:n_orb, :n_orb]

# Get Delta(iw) and chemical_potential from G0(iw)
def get_h0_Delta(G0_iw):
    h0_lst, Delta_iw = [], G0_iw.copy()
    for bl in G0_iw.indices:
        Delta_iw[bl] << iOmega_n - inverse(G0_iw[bl])
        tail, err = fit_hermitian_tail(Delta_iw[bl])
        Delta_iw[bl] << Delta_iw[bl] - tail[0]
        h0_lst.append(tail[0])
    return h0_lst, Delta_iw

h0_lst, Delta_iw = get_h0_Delta(G0_iw)
chemical_potential = []
for h0 in h0_lst:
    chemical_potential += [-l for l in linalg.eig(h0)[0].real]

# Fourier-transform to get Delta(tau)
tau_mesh = MeshImTime(beta, 'Fermion', n_tau)
Delta_tau = BlockGf(mesh = tau_mesh, gf_struct = gf_struct)
for name, g0 in Delta_iw:
    known_moments = make_zero_tail(Delta_iw[name], 1)
    tail, err = fit_hermitian_tail(Delta_iw[name], known_moments)
    Delta_tau[name] << make_gf_from_fourier(Delta_iw[name], tau_mesh, tail).real
###########################################################

# --------- Construct the CTSEG solver ----------
constr_params = {
    'beta': beta,
    'gf_struct': gf_struct,
    'n_tau': n_tau,
}
S = Solver(**constr_params)

# Input Delta(tau)
S.Delta_tau << Delta_tau

# Impurity interaction Hamiltonian
Umat, Upmat = U_matrix_kanamori(n_orb, U_int=U, J_hund=J)
h_int = h_int_density(block_names, n_orb, Umat, Upmat, off_diag=True)

# --------- Solve parameters ----------
solve_params = {
    'h_int': h_int,
    'chemical_potential': chemical_potential,
    'n_warmup_cycles': 5000,
    'n_cycles': 50000,
    'length_cycle': 100,
}

# ---------- Solve ----------
S.solve(**solve_params)

# --------- Save output ---------
if mpi.is_master_node():
    with h5.HDFArchive("dimer.out.h5", 'w') as A:
        A['G_tau'] = S.results.G_tau

# --------- Compare to reference ----------      
    #h5diff("dimer.out.h5", "dimer.ref.h5", precision=1e-8)
    

