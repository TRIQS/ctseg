# Copyright (c) 2022-2024 Simons Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You may obtain a copy of the License at
#     https:#www.gnu.org/licenses/gpl-3.0.txt
#
# Authors: Olivier Parcollet, nkavokine

from .solver_core import SolverCore

from triqs.gf import *
from triqs.utility import mpi

# === The SolverCore Wrapper

class Solver(SolverCore):
    """
    The solver class.
    """

    def __init__(self, **kwargs):
        """
        Initialize the solver.

        Parameters
        ----------
        .. include:: ../../python/triqs_ctseg/parameters_constr_params_t.rst
        """

        kwargs['gf_struct'] = fix_gf_struct_type(kwargs['gf_struct'])

        # Initialize the solver
        SolverCore.__init__(self, **kwargs)

    def solve(self, **kwargs):
        """
        Solve the impurity problem.

        Parameters
        ----------
        .. include:: ../../python/triqs_ctseg/parameters_solve_params_t.rst
        """

        # Solve the impurity problem
        solve_status = SolverCore.solve(self, **kwargs)
        return solve_status
