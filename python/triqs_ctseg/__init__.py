# Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE.txt in the root of this distribution for details.

r"""
The CTSEG python module. 
"""
from .solver import Solver
from .solver_core import SolverCore

__all__ = ['Solver', 'SolverCore']


class Cpp2pyInfo:
    table_imports = {}
    table_converters = {}
