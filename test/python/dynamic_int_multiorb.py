import numpy as np
from h5 import HDFArchive

from triqs.gf import make_gf_imtime
from triqs.operators import util

from triqs_ctseg import SolverCore as Solver

# load input data from h5
# the input data is generated from the github.com/TRIQS/solid_dmft/tree/3.3.x/test/python/svo_gw_emb_dyn test
with HDFArchive('dynamic_int_multiorb_input.h5', 'r') as ar:
    delta_tau = ar['delta_tau']
    chemical_potential = ar['chemical_potential']
    Uloc_dlr_2idx = ar['Uloc_dlr_2idx']
    Uloc_dlr_2idx_prime = ar['Uloc_dlr_2idx_prime']
    Vloc = ar['Vloc']

# Number of orbitals and time mesh for interaction
n_orb = Vloc.shape[0]
n_tau_k = 10001

# Solver construction parameters
constr_params = {
    'gf_struct': [(block, 1) for block in delta_tau.indices],
    'beta': delta_tau.mesh.beta,
    'n_tau': len(delta_tau.mesh),
    'n_tau_k': n_tau_k,
}

# Construct solver
S = Solver(**constr_params)

# Hybridization Delta(tau)
for name, block in S.Delta_tau:
    block << delta_tau[name]

# Interaction Hamiltonian
Umat, Upmat = util.reduce_4index_to_2index(Vloc)

h_int = util.h_int_density(['down', 'up'], n_orb, off_diag=False, U=Umat, Uprime=Upmat)

# create full frequency interaction on full time mesh
Uloc_tau_2idx = make_gf_imtime(Uloc_dlr_2idx, n_tau=n_tau_k)
Uloc_tau_2idx_prime = make_gf_imtime(Uloc_dlr_2idx_prime, n_tau=n_tau_k)

# fill D0_tau from Uloc_tau_2idx and Uloc_tau_2idx_prime
# same spin interaction
S.D0_tau[0:n_orb, 0:n_orb] << Uloc_tau_2idx.real
S.D0_tau[n_orb:2*n_orb, n_orb:2*n_orb] << Uloc_tau_2idx.real
# opposite spin interaction
S.D0_tau[0:n_orb, n_orb:2*n_orb] << Uloc_tau_2idx_prime.real
S.D0_tau[n_orb:2*n_orb, 0:n_orb] << Uloc_tau_2idx_prime.real

# Solve parameters
solve_params = {
    'h_int': h_int,
    'chemical_potential': chemical_potential,
    'length_cycle': 600,
    'n_warmup_cycles': 2000,
    'n_cycles': 40000,
    'measure_perturbation_order_histograms': True,
    'measure_statehist': True,
    'measure_n': True
}

# Solve
S.solve(**solve_params)

# very basic check of the results
ref_occ = np.array([0.1673 for occ in range(2*n_orb)])
assert(np.allclose(S.results.densities, ref_occ, atol=1e-2))

