from triqs.gf import *
from triqs_ctseg import SolverCore as Solver
import triqs.utility.mpi as mpi
from triqs.gf.descriptors import Function
from triqs.plot.mpl_interface import *
import numpy as np
import h5

with h5.HDFArchive("/mnt/home/nkavokine/ctseg_J/build/test/c++/anderson.ref.h5", 'r') as Af:
    gref = Af["(ctqmc.G_tau()[0])"]

oplot(gref)

beta = 20
gf_struct = [("up", [0]), ("down", [0])]
S = Solver(beta=beta,  gf_struct=gf_struct)
