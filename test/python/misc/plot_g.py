from triqs.gf import *
from triqs.plot.mpl_interface import oplot
from h5 import *
A = HDFArchive("single_bath_seg.output.h5")
a=A["G_tau"]["up"].data
import numpy
numpy.savetxt('gtau.dat',a[:,0,0])
a=A["K_tau"].data
import numpy
numpy.savetxt('ktau.dat',a[:,0,0])
