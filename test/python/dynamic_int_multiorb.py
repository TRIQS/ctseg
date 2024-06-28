import numpy as np
from h5 import HDFArchive

from triqs.gf import make_gf_imtime
from triqs.operators import util
from triqs.operators import n
from triqs.utility import mpi
from triqs.utility.h5diff import h5diff


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
n_tau_bosonic = 10001

# Solver construction parameters
gf_struct = [(block, 1) for block in delta_tau.indices]
constr_params = {
    'gf_struct': gf_struct,
    'beta': delta_tau.mesh.beta,
    'n_tau': len(delta_tau.mesh),
    'n_tau_bosonic': n_tau_bosonic,
}
print([(block, 1) for block in delta_tau.indices])

# Construct solver
S = Solver(**constr_params)

# Hybridization Delta(tau)
for name, block in S.Delta_tau:
    block << delta_tau[name]

# Interaction Hamiltonian
Umat, Upmat = util.reduce_4index_to_2index(Vloc)

h_int = util.h_int_density(['down', 'up'], n_orb, off_diag=False, U=Umat, Uprime=Upmat)

# Quadratic Local Hamiltonian
h_loc0 = sum(-chemical_potential[i] * n(block, 0) for i, (block, s) in enumerate(gf_struct))


# create full frequency interaction on full time mesh
Uloc_tau_2idx = make_gf_imtime(Uloc_dlr_2idx, n_tau=n_tau_bosonic)
print(Uloc_tau_2idx)
Uloc_tau_2idx_prime = make_gf_imtime(Uloc_dlr_2idx_prime, n_tau=n_tau_bosonic)

# fill D0_tau from Uloc_tau_2idx and Uloc_tau_2idx_prime
# same spin interaction
for a in range(n_orb):
    for b in range(n_orb):
        S.D0_tau["up_{}".format(a), "up_{}".format(b)][0, 0] << Uloc_tau_2idx[a, b].real
        S.D0_tau["down_{}".format(a), "down_{}".format(b)][0, 0] << Uloc_tau_2idx[a, b].real
# opposite spin interaction
for a in range(n_orb):
    for b in range(n_orb):
        S.D0_tau["up_{}".format(a), "down_{}".format(b)][0, 0] << Uloc_tau_2idx_prime[a, b].real
        S.D0_tau["down_{}".format(a), "up_{}".format(b)][0, 0] << Uloc_tau_2idx_prime[a, b].real

# Solve parameters
solve_params = {
    'h_int': h_int,
    'h_loc0': h_loc0,
    'length_cycle': 100,
    'n_warmup_cycles': 2000,
    'n_cycles': 10000,
    'measure_pert_order': True,
    'measure_state_hist': True,
    'measure_densities': True,
    'measure_F_tau': True,
    'measure_nn_static': True,
    'measure_nn_tau': True,

}

# Solve
S.solve(**solve_params)

# Save and compare to reference
if mpi.is_master_node():
    with HDFArchive("dynamic_int_multiorb.out.h5", 'w') as A:
        A['G_tau'] = S.results.G_tau
        A['F_tau'] = S.results.F_tau
        A["nn_tau"] = S.results.nn_tau
        A['nn'] = S.results.nn_static
        A['densities'] = S.results.densities

    h5diff("dynamic_int_multiorb.out.h5", "dynamic_int_multiorb.ref.h5", precision=1e-9)

