import sys
import importlib
import numpy as np
import matplotlib.pyplot as plt
from triqs.gf import *
import triqs.utility.mpi as mpi
from triqs.gf.descriptors import Function
from triqs.operators import *
from triqs.plot.mpl_interface import oplot
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

beta = 20
U = 2.0
#mu   = 0.2
mu = U/2
epsilon = 0.3
l = 1.0  # electron boson coupling
w0 = 1.0  # screening frequency
n_tau = 1001
n_iw = 200
random_seed = 23488
length_cycle = 1
n_warmup_cycles = 1
n_cycles = 1
measure_gt = True
measure_n = True
move_move = False

has_Dt = True

# Local Hamiltonian
H = U*n("up", 0)*n("down", 0)

# Hybridization
delta = GfImFreq(indices=[0], beta=beta, n_points=n_iw)
delta << inverse(iOmega_n - epsilon)  # + inverse(iOmega_n + epsilon)

# Parameters
p = {}
p["random_seed"] = random_seed
p["length_cycle"] = length_cycle
p["n_warmup_cycles"] = n_warmup_cycles
p["n_cycles"] = n_cycles
p["measure_gt"] = measure_gt
p["measure_n"] = measure_n
p["move_insert_segment"] = True
p["move_remove_segment"] = True
p["move_split_segment"] = True
p["move_regroup_segment"] = True
p["move_move_segment"] = move_move
p["hartree_shift"] = [mu, mu]

# Solver!
Snew = ctseg_new.SolverCore(beta=beta,
                            gf_struct=[['down', 1], ['up', 1]],
                            n_tau=n_tau,
                            n_tau_k=n_tau,
                            )


Snew.Delta_tau << Fourier(delta)

if has_Dt:
    D0w = GfImFreq(indices=[0], beta=beta, statistic='Boson', n_points=n_iw)
    D0w << Function(lambda w: l**2*2*w0/(w**2-w0**2))
    Snew.D0_tau << Fourier(D0w)

Snew.solve(h_int=H, **p)
