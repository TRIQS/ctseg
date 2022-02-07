from triqs_ctseg.scripts.subohmic import *
from triqs_ctseg.block_matrix import *

beta=10.0
Vsq=0.1
g=0.1
l=0.1
time_max=-1
n_cycles=5000
filename="spinspin_subohmic"
n_w=50

run_spinspin(beta, Vsq, g, time_max, rank=-1, from_scratch=True, U=0.4, l=l, measure_fw=False, split_regroup_v1=False,split_regroup_v2=False, sign=1.0, mu=0.2, move_move=False, move_swap_empty_lines=True, measure_gw=False, measure_g2w=False, n_w=n_w, filename=filename, n_cycles=n_cycles)

from triqs.utility.h5diff import h5diff
h5diff(filename+".out.h5", filename+".ref.h5")
