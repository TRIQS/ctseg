# Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE.txt in the root of this distribution for details.

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
