# Copyright (c) 2022--present, The Simons Foundation
# Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

r"""User-facing CTSEG solver.

This module exposes :class:`Solver`, a thin Python wrapper around 
:class:`~triqs_ctseg.solver_core.SolverCore`. 
"""

from .solver_core import SolverCore, ConstrParamsT, SolveParamsT

from triqs.gfs import *
from triqs.utility import mpi

# === The SolverCore Wrapper

class Solver(SolverCore):
    r"""Continuous-time hybridization-expansion impurity solver in the segment picture.

    Thin Python wrapper around :class:`~triqs_ctseg.solver_core.SolverCore`. 
    After construction the user is expected to assign the hybridisation function 
    :attr:`Delta_tau` and, optionally, the retarded density-density interaction 
    :attr:`D0_tau` and the transverse spin-spin interaction :attr:`Jperp_tau` before 
    calling :meth:`solve`. 
    
    Measured observables are read from the inherited :attr:`results` member (see 
    :class:`~triqs_ctseg.solver_core.ResultsT`).

    Parameters
    ----------
    **kwargs
        Construction parameters forwarded to
        :class:`~triqs_ctseg.solver_core.ConstrParamsT`; see that class for
        the full list with defaults.
    """

    def __init__(self, **kwargs):
        r"""Initialise the solver.

        Parameters
        ----------
        **kwargs
            Construction parameters forwarded to
            :class:`~triqs_ctseg.solver_core.ConstrParamsT`; see that class
            for the full list with defaults.
        """

        kwargs['gf_struct'] = fix_gf_struct_type(kwargs['gf_struct'])

        # Initialize the solver
        SolverCore.__init__(self, ConstrParamsT(**kwargs))

    def solve(self, **kwargs):
        r"""Solve the impurity problem.

        Runs the Monte Carlo simulation with the requested moves and measures
        and stores the accumulated observables in :attr:`results`.

        Parameters
        ----------
        **kwargs
            Solve parameters forwarded to
            :class:`~triqs_ctseg.solver_core.SolveParamsT`; see that class
            for the full list with defaults.
        """

        # Solve the impurity problem
        solve_status = SolverCore.solve(self, SolveParamsT(**kwargs))
        return solve_status

    @property
    def density_matrix(self):
        r"""TTI diagonal impurity density matrix in the occupation basis.

        The value is stored in the legacy ``results.state_hist`` HDF field for
        backward compatibility.
        """
        return self.results.state_hist
