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

r"""
The CTSEG python module. 
"""
from .solver import Solver
from .solver_core import SolverCore

__all__ = ['Solver', 'SolverCore']


class Cpp2pyInfo:
    table_imports = {}
    table_converters = {}
