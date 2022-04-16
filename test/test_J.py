# %%
import sys
import importlib
import numpy as np
from triqs.gf import *
import triqs.utility.mpi as mpi
from triqs.gf.descriptors import Function
from triqs.operators import *
import h5

NEW_PATH = "/mnt/home/nkavokine/ctseg_J/build/python/triqs_ctseg/__init__.py"
OLD_PATH = "/mnt/home/nkavokine/ctseg_OLD/build/python/triqs_ctseg/__init__.py"

spec_new = importlib.util.spec_from_file_location("triqs_ctseg_J", NEW_PATH)
ctseg_new = importlib.util.module_from_spec(spec_new)
sys.modules[spec_new.name] = ctseg_new
spec_new.loader.exec_module(ctseg_new)
spec = importlib.util.spec_from_file_location("triqs_ctseg", OLD_PATH)
ctseg_old = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ctseg_old
spec.loader.exec_module(ctseg_old)


# %% [markdown]
# Parameters

# %%
beta = 10
U = 4.0
mu = U/2
J = 0.5
n_tau = 2051
n_tau_k = 2001
n_iw = 1025
#random_seed = 34788 + 928374
length_cycle = 50
n_warmup_cycles = 10000
n_cycles = 1000000  # 7500000
measure_gt = True
measure_ft = True
measure_n = True
measure_nnt = False
move_move = True
dynamical_U = True
spin = True
old_code = False
new_code = True

with h5.HDFArchive("/mnt/home/nkavokine/pyscripts/solver-jperp-ctseg-ctint-test/run-ctseg-Jperp1-J0.50-beta10.00-U4.00-muS1.00.out.h5", 'r') as Af:
    g0 = Af["dmft_loop/i_001/S/G0_iw/up"]
    Q_tau = Af["dmft_loop/i_000/Q_tau"]

# Local Hamiltonian
H = U*n("up", 0)*n("down", 0)

# Hybridization
delta = GfImFreq(indices=[0], beta=beta, n_points=n_iw)
invg0 = GfImFreq(indices=[0], beta=beta, n_points=n_iw)
invg0 << inverse(g0)
delta << iOmega_n + mu - invg0

# Hybridization
#epsilon = 0.2
#delta = GfImFreq(indices=[0], beta=beta, n_points=n_iw)
# delta << inverse(iOmega_n - epsilon)  # + inverse(iOmega_n + epsilon)

if new_code:

    # %% [markdown]
    # Setup new solver

    # %%
    # Parameters
    p = {}
    #p["random_seed"] = random_seed
    p["length_cycle"] = length_cycle
    p["n_warmup_cycles"] = n_warmup_cycles
    p["n_cycles"] = n_cycles
    p["measure_gt"] = measure_gt
    p["measure_ft"] = measure_ft
    p["measure_n"] = measure_n
    p["measure_nnt"] = measure_nnt
    p["hartree_shift"] = [mu, mu]

    p["move_insert_segment"] = True
    p["move_remove_segment"] = True
    p["move_split_segment"] = True
    p["move_regroup_segment"] = True
    p["move_insert_spin_segment"] = True
    p["move_remove_spin_segment"] = True
    p["move_split_spin_segment"] = True
    p["move_regroup_spin_segment"] = True
    p["move_swap_spin_lines"] = True
    p["move_move_segment"] = move_move

    # Solver!
    Snew = ctseg_new.SolverCore(beta=beta,
                                gf_struct=[['down', 1], ['up', 1]],
                                n_tau=n_tau,
                                n_tau_k=n_tau_k,
                                )

    Snew.Delta_tau << Fourier(delta)

    D0 = GfImTime(indices=[0, 1], beta=beta,
                  statistic='Boson', n_points=n_tau_k)
    D0[0, 0] = -Q_tau[0, 0]
    D0[1, 1] = -Q_tau[0, 0]
    D0[0, 1] = Q_tau[0, 0]
    D0[1, 0] = Q_tau[0, 0]
    if spin:
        Snew.Jperp_tau << -J**2*Q_tau
    if dynamical_U:
        Snew.D0_tau << 0.25*J**2*D0
    # %%
    Snew.solve(h_int=H, **p)
    if mpi.is_master_node():
        with h5.HDFArchive("test_J_new.h5", 'w') as A:
            A['G_tau'] = Snew.results.G_tau
            A['F_tau'] = Snew.results.F_tau


if old_code:
    # %% [markdown]
    # Setup old solver

    # %%
    # Parameters
    p = {}
    #p["random_seed"] = random_seed
    p["length_cycle"] = length_cycle
    p["n_warmup_cycles"] = n_warmup_cycles
    p["n_cycles"] = n_cycles
    p["measure_gt"] = measure_gt
    p["measure_ft"] = measure_ft
    p["measure_nnt"] = measure_nnt
    p["move_insert_segment"] = True
    p["move_remove_segment"] = True
    p["move_move"] = False
    p["hartree_shift"] = [0, 0]

    gf_struct = [("up", 1), ("down", 1)]

    # Solver!
    Sold = ctseg_old.SolverCore(beta=beta,
                                gf_struct=gf_struct,
                                n_tau=n_tau,
                                n_tau_k=n_tau,
                                n_tau_jperp=n_tau,
                                n_iw=n_iw
                                )

    Sold.G0_iw << g0
    #Sold.G0_iw << inverse(iOmega_n + mu - delta)
    if dynamical_U:
        Sold.D0_iw['up|up'] << Fourier(-0.25*J**2*Q_tau)
        Sold.D0_iw['down|down'] << Fourier(-0.25*J**2*Q_tau)
        Sold.D0_iw['up|down'] << Fourier(0.25*J**2*Q_tau)
        Sold.D0_iw['down|up'] << Fourier(0.25*J**2*Q_tau)
    if spin:
        Sold.Jperp_iw << Fourier(-J**2*Q_tau)

    # %%
    Sold.solve(h_int=H, **p)
    if mpi.is_master_node():
        with h5.HDFArchive("test_J_old.h5", 'w') as A:
            A['G_tau'] = Sold.G_tau
            A['F_tau'] = Sold.F_tau
