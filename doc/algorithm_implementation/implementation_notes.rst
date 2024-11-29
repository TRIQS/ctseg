.. _implementation_notes:

Implementation notes 
====================

*A guide to reading the code.*

CTSEG is an application based on the `app4triqs <https://triqs.github.io/app4triqs/unstable/index.html>`_ skeleton. The code is written in C++, and a python wrapper
is generated automatically using ``cpp2py``. The user interacts with the solver in python. 

Solver core
***********

The solver is an instance of the ``solver_core`` class. It is constructed based on parameters (``constr_params``) passed as
a python dictionary (see ``params.hpp``). It contains the set of inputs (``inputs.hpp``), 
the set of results (``results.hpp``) and the ``solve`` method. The latter runs 
the QMC calculation for a given set of inputs, and given the additional ``solve_params`` (see ``params.hpp``). 

The solve method works with three structures: ``params_t p`` (the set of all parameters), ``configuration_t config``
(the Monte Carlo configuration, see below) and ``work_data_t wdata``, that contains all the data and auxiliary 
methods that are used by the Monte Carlo moves and measurements. It further defines an instance of the 
`mc_generic <https://triqs.github.io/triqs/latest/documentation/manual/triqs/mc_tools/contents.html>`_ class 
(a component of TRIQS) that takes care of implementing the Metropolis algorithm. The 
``mc_generic CTQMC`` object is supplied with "moves" and "measures": classes that satisfy a particular concept. Its 
``warmup_and_accumulate`` method then takes care of choosing moves at random, accepting or rejecting them according 
to the Metropolis criterion, and sampling the observables every few moves. When the specified number of Monte Carlo 
cycles has been completed, ``collect_results`` averages the results across nodes and makes them accessible to the user
through the ``results`` structure. 

Inputs
******

The inputs are the numerical data required for defining the impurity problem. They are stored in the ``inputs``
structure, but they can be accessed directly thanks to dedicated methods of the ``solver_core`` class. 

* ``Delta_tau`` is the hybridization function. It should be supplied as a TRIQS ``BlockGf``. 
  The choice of the block structure is usually not unique, but using as many blocks as possible ensures the fastest 
  determinant updates. 
* ``D0_tau`` is the dynamical part of the density-density interaction. Note that the static part is supplied in 
  the ``solve_params`` as part of the interaction Hamiltonian. It should be supplied as a TRIQS ``Block2Gf``
  (since it is no necessarily block-diagonal like ``Delta_tau``).
* ``Jperp_tau`` is the perpendicular spin-spin interaction. It should be supplied as a :math:`1 \times 1` matrix-valued
  ``GfImTime``.


Configuration
*************

The basic units of a configuration (see ``configuration.hpp``) are segments and :math:`J_{\perp}` lines. A segment 
is a structure containing the time corresponding to a :math:`c` operator, the time corresponding to a :math:`c^{\dagger}`
operator, and a boolean for each operator indicating whether it is attached to a :math:`J_{\perp}` line. 

.. note::

    By convention, we orient the imaginary time axis from :math:`\beta` on the left to 0 on the right. 

A segment may be cyclic if :math:`\tau_c < \tau_{c^{\dagger}}`. A line that is occupied at all times
(a zeroth order term in the hybridization expansion) is represented by a segment with 
:math:`\tau_{c^{\dagger}} = 0` and :math:`\tau_c = \beta`. 

A :math:`J_{\perp}` line is a structure containing the time corresponding to an :math:`S^+` operator
and the time corresponding to an :math:`S^-` operator. No orbital indices are stored because the 
:math:`J_{\perp}` expansion is only implemented for a single orbital. 

The structure of the configuration is inherited from the structure of the hybridization function. The 
hybridization function is matrix-valued, and its line (or column) indices are termed *colors*. The configuration
consists of a list (``std::vector``) of lists of segments (one per color), and a list of :math:`J_{\perp}` lines. 
The block structure of the hybridization function is irrelevant for the configuration (it is only used for the determinants, see below). 

.. warning::

    The segment lists are sorted so that the :math:`c` operators are in **decreasing** time order. There is 
    no imposed ordering for the list of :math:`J_{\perp}` lines. When a :math:`J_{\perp}` expansion is carried 
    out, it is assumed that color 0 is spin up and color 1 is spin down. 

Imaginary time points
*********************

Imaginary time points are represented as 64 bit integers, via the custom type ``tau_t`` (see ``tau_t.hpp``). This discretization
on a very fine grid allows for exact comparisons, which are notoriously dangerous with floating point numbers. 
The Monte Carlo moves explicitly forbid degenerate segments with the :math:`c` and :math:`c^{\dagger}` operators at 
equal times, as well as coinciding operators for segments of the same color. The ``tau_t`` type also implements 
:math:`\beta`-periodicity for the addition of imaginary times. 

Work data
*********

The ``work_data`` structure (see ``work_data.hpp``) contains data and methods that are used by the Monte Carlo moves 
and measurements. Most importantly, its construction involves computing the dynamical interaction kernel :math:`K(\tau)`, 
and initializing the determinant for every block of the hybridization matrix :math:`[\Delta]`. ``work_data.hpp`` also 
contains auxiliary functions for the Monte Carlo moves: in particular, ``trace_sign``, that computes the sign of the trace 
from the times of the hybridized operators stored in the ``dets`` object. 

Determinants
************

The hybridization matrix :math:`[\Delta]` is stored in ``work_data`` as a list (``std::vector``) of its blocks. 
Each block is a ``det_manip`` object. The ``det_manip`` class is a part of TRIQS, which implements
a fast computation of the change to :math:`\mathrm{det}[\Delta]` upon insertion or removal of a line 
and column in :math:`[\Delta]`. The ``det_manip`` object is constructed from the input hybridization function 
via an "adaptor" that specifies how it is to be called (see ``dets.hpp``). The alias ``det_t`` is used for the 
type of :math:`[\Delta]`. 

.. warning::

    Within a given block, the lines and columns of the hybridization matrix :math:`[\Delta]` are arranged in **increasing**
    order of the corresponding times, regardless of color. 

Miscellaneous
*************

* ``util.hpp`` contains utility functions used by the Monte Carlo moves and measurements. 
* ``logs.hpp`` contains the preprocessor directives that set the logging level at compile time, and macros that facilitate logging 
  (we use the ``spdlog`` library). 
* ``invariants.hpp`` contains functions that check invariants of the algorithm 
  at every move if the code is compiled in debug mode. Checks include the time-ordering of the segments and of the lines and 
  columns of :math:`[\Delta]`, the consistency between the times stored in the configuration and in :math:`[\Delta]`, and the 
  consistency between the :math:`J_{\perp}` lines and the labels stored in the segments. 

Automatic tests
***************

Several automatic tests are supplied with CTSEG. They run short simulations to assess whether the software functions correctly.
* ``anderson.py``, ``dynamical_U.py``, ``Jperp.py``, ``spin_spin.py`` and ``multiorb.py`` are non-regression tests, that compare their results to a reference generated by CTSEG itself. 
  They cover most situations that might be handled by CTSEG: single or multiple orbitals, static, dynamic or spin-spin interactions. 
  In addition, two tests written directly in C++, ``anderson.cpp`` and ``spin_spin.cpp``, are supplied for developers who may wish to check for non-regression without re-generating the Python wrapper. 
* ``dimer.py`` may be used to check the correctness of CTSEG on an exactly solvable problem (two impurity sites coupled to two bath sites). The reference data can be found
  in ``dimer_pyed.ref.h5``. A meaningful comparison to the reference requires a long solver run (for this particular problem, the result can have significant fluctuations between machines for short runs). 
  Such a long run is not carried out as part of the automatic testing procedure and should be set up by the user if an in-depth assessment of correctness is required. 
* Additionally, the correctness of CTSEG for problems involving dynamical and spin-spin interaction can be checked by comparing the output of ``spin_spin.py`` to the data in ``ctint.ref.h5``, 
  generated by a very long CTINT run. CTINT carries out an altogether different expansion to solve the problem: therefore, agreement between CTSEG and CTINT is a good indication for the correctness 
  of both solvers. 
