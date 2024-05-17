r"""
The CTSEG python module. 
"""
from .solver import Solver
from .solver_core import SolverCore

__all__ = ['Solver', 'SolverCore']


class Cpp2pyInfo:
    table_imports = {}
    table_converters = {}
