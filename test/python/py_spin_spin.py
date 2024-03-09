# Single orbital with dynamical spin-spin interactions. 
# Can be compared with reference obtained with CTINT (ctint.ref.h5)
from triqs.gf import *
import triqs.utility.mpi as mpi
from triqs.gf.descriptors import Function
from triqs.operators import *
import h5
from triqs.utility.h5diff import h5diff
from triqs_ctseg import SolverCore as Solver

# Parameters
beta = 10
U = 4.0
mu = U/2
J = 0.5

# Local Hamiltonian
H = U*n("up", 0)*n("down", 0)

# Get inputs from reference file
with h5.HDFArchive("ctint.ref.h5", 'r') as Af:
    g0 = Af["dmft_loop/i_001/S/G0_iw/up"]
    Q_tau = Af["dmft_loop/i_000/Q_tau"]

# Hybridization
n_iw = 1025
delta = GfImFreq(indices=[0], beta=beta, n_points=n_iw)
invg0 = GfImFreq(indices=[0], beta=beta, n_points=n_iw)
invg0 << inverse(g0)
delta << iOmega_n + mu - invg0

# Parameters
p = {}
p["length_cycle"] = 50
p["n_warmup_cycles"] = 1000
p["n_cycles"] = 10000
p["measure_nnt"] = True
p["measure_nn"] = True
p["measure_ft"] = True
p["hartree_shift"] = [mu, mu]

# Construct solver
n_tau = 2051
n_tau_k = 2001
S = Solver(beta=beta,
            gf_struct=[('down', 1), ('up', 1)],
            n_tau=n_tau,
            n_tau_k=n_tau_k,
            )

# Solver inputs
S.Delta_tau << Fourier(delta)

D0 = GfImTime(indices=[0, 1], beta=beta, statistic='Boson', n_points=n_tau_k)
D0[0, 0] = -Q_tau[0, 0]
D0[1, 1] = -Q_tau[0, 0]
D0[0, 1] = Q_tau[0, 0]
D0[1, 0] = Q_tau[0, 0]
S.Jperp_tau << -J**2*Q_tau
S.D0_tau << 0.25*J**2*D0

# Solve parameters
p = {}
p["length_cycle"] = 50
p["n_warmup_cycles"] = 1000
p["n_cycles"] = 10000
p["measure_nnt"] = True
p["measure_nn"] = True
p["measure_ft"] = True
p["hartree_shift"] = [mu, mu]

# Solve!
S.solve(h_int=H, **p)

if mpi.is_master_node():
    with h5.HDFArchive("py_spin_spin.out.h5", 'w') as A:
        A['G_tau'] = S.results.G_tau
        A['F_tau'] = S.results.F_tau
        A['nn_tau'] = S.results.nn_tau
        A['nn'] = S.results.nn_static
        A['densities'] = S.results.densities

    h5diff("py_spin_spin.out.h5", "py_spin_spin.ref.h5", precision=1e-9) 
