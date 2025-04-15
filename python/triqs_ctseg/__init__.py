# Copyright (c) 2022--present, The Simons Foundation
# Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

r"""
The CTSEG python module. 
"""
from .solver import Solver
from .solver_core import SolverCore

__all__ = ['Solver', 'SolverCore']


class Cpp2pyInfo:
    table_imports = {}
    table_converters = {}
