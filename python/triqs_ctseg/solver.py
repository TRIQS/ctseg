from .solver_core import SolverCore

from triqs.gf import *
from triqs.utility import mpi

# === The SolverCore Wrapper

class Solver(SolverCore):

    def __init__(self, **kwargs):
        """
        Initialise the solver.

        Parameters
        ----------
        Cf. C++ documentation of SolverCore
        """

        kwargs['gf_struct'] = fix_gf_struct_type(kwargs['gf_struct'])

        # Initialise the core solver
        SolverCore.__init__(self, **kwargs)

    def solve(self, **kwargs):
        """
        Solve the impurity problem.

        Parameters
        ----------
        Cf. C++ documentation of SolverCore.solve(..)
        """

        # Solve the impurity problem
        solve_status = SolverCore.solve(self, **kwargs)
        return solve_status
