.. _documentation:

Guided tour
=========================

The 'cthyb-segment' solver is a continuous-time Quantum Monte Carlo solver aimed at solving an :doc:`Anderson impurity model <definitions>` with static and dynamical density-density interactions.

Below is a very short example demonstrating the basic usage of the solver.

.. runblock:: python

  from triqs.gf import *
  from h5 import *
  from triqs.operators import *
  import triqs.utility.mpi as mpi
  from triqs_ctseg import Solver

  #Green's function: two 1x1 blocks
  gf_struct = [("up",[0]), ("down",[0])]
  #initialize the solver
  S = Solver(
       beta = 10.0, #inverse temperature
       gf_struct = gf_struct )

  #set up the static interaction
  U = 2.0 #static interaction

  #impurity non-interacting Green's function
  mu = 1.0 #chemical potential
  eps0 = 0.0 #fermionic bath level
  V = 0.5 #coupling to fermionic level
  S.G0_iw << \
     inverse(iOmega_n + mu \
             -V**2 * inverse(iOmega_n - eps0)\
             -V**2 * inverse(iOmega_n + eps0))

  #dynamical charge interactions
  w0 = 1.0 #bosonic bath level
  lambd = 0.4 #coupling to bosonic level
  from itertools import product
  for b1,b2 in product(dict(gf_struct).keys(), repeat=2):
   S.D0_iw[b1+"|"+b2] << \
     lambd**2*(inverse(iOmega_n - w0) \
                   - inverse(iOmega_n + w0))

  #run the Monte-Carlo simulation
  S.solve(
   h_int = U * n('up',0)*n('down',0), #local Hamiltonian
   n_cycles  = 10000, #Monte-Carlo cycles
   length_cycle = 100, #length of a cycle
   n_warmup_cycles = 1000,#thermalization cycles
   n_w_b_vertex = 2, #bosonic freqs. in vertex
   n_w_f_vertex = 20,#fermionic freqs. in vertex
   measure_nn = True, #density-density corr.
   measure_nnw = True, #density-density corr. fct
   measure_gw = True, #one-particle Green fct
   measure_fw = True, #one-particle improved est.
   measure_g2w = True,#3-leg corr.function
   measure_g3w = True,#4-leg corr. function
  )

  #store results in a HDF5 archive
  if mpi.is_master_node():
   with HDFArchive("ct_seg.h5",'a') as A:
    A['G_tau'] = S.G_tau
    A['nn'] = S.nn
    A['nn_iw'] = S.nn_iw
    A['Sigma_iw']  = S.Sigma_iw
    A["G_2w"] = S.G_2w
    A["G_3w"] = S.G_3w



To use the solver, one has to go through the following steps:

*Import the* ``Solver`` *class*::

  from triqs_ctseg import Solver 
 
*Declare and construct a* ``Solver`` *instance*::

  gf_struct = [("up",[0]), ("down",[0])]
  S = Solver(
       beta = 10.0, #inverse temperature
       gf_struct = gf_struct )

Here, ``gf_struct`` is a Python dictionary that specifies the block structure of the Green's function, here one :math:`1\times 1` (``[0]`` is of length 1) for each spin (``up``, ``down``).

Some construction parameters are mandatory, others are optional. See :doc:`here <solver_core>` for a complete list.
 
*Initialize the inputs of the solver*. In this script, we want to define the noninteracting Green's function as:

.. math::

  \begin{equation}
  \left[\mathcal{G}_{0}\right]^{\sigma}(i\omega_{n})=\left[i\omega_{n}-\varepsilon_{a\sigma}-\Delta_{ab}^{\sigma}(i\omega_{n})\right]^{-1},\label{eq:G0_def}
  \end{equation}

with 

.. math::

  \begin{align*}
  \Delta^\sigma(i\omega_{n}) & =\frac{V^{2}}{i\omega_{n}-\varepsilon_{0}} +\frac{V^{2}}{i\omega_{n}+\varepsilon_{0}},\\
  \end{align*}

In the script, this is done via an *accessor* called ``G0_iw``. It is set by the user in the following way::

  S.G0_iw << \
     inverse(iOmega_n + mu \
             -V**2 * inverse(iOmega_n - eps0)\
             -V**2 * inverse(iOmega_n + eps0))

Similarly, we want to define dynamical interactions:

.. math::

  \begin{align*}
  \mathcal{D}_{0}^{\sigma,\sigma'}(i\Omega_{m}) & =\frac{2\lambda^{2}\omega_0}{(i\Omega_{m})^{2}-\omega_{0}^2}.
  \end{align*}


This is done through the ``D0_iw`` accessor::

  for b1,b2 in product(dict(gf_struct).keys(), repeat=2):
   S.D0_iw[b1+"|"+b2] << \
     lambd**2*(inverse(iOmega_n - w0) \
                   - inverse(iOmega_n + w0))
 
Other accessors are used to read out the observables after running the solver, like at the end of script::

    A["G_tau"] = S.G_tau
    A['nn'] = S.nn
    A['nn_iw'] = S.nn_iw
    A['Sigma_iw']  = S.Sigma_iw
    A["G_2w"] = S.G_2w
    A["G_3w"] = S.G_3w

Here, ``G_tau`` is the one-particle Green's function :math:`G^\sigma(\tau)`, ``nn`` the density-density correlator :math:`\chi^{\sigma\sigma'}_{ab}`, ``nn_iw`` the density-density correlation function :math:`\chi^{\sigma\sigma'}_{ab}(i\Omega)`, ``Sigma_iw`` the self-energy :math:`\Sigma^\sigma(i\omega)`, ``G_2w`` the three-point correlation function :math:`\chi^{\sigma\sigma'}_{abc}(i\omega,i\Omega)` and ``G_3w`` the four-point correlation function :math:`G^{2,\sigma\sigma'}_{abcd}(i\omega,i\omega',i\Omega)`. All these observables are defined :doc:`here <definitions>` or :doc:`there <measures>`.

For a complete list of the accessors (both input and output), see :doc:`here <solver_core>`.

*Call the* ``solve()`` *method*::

  S.solve(
   h_int = U * n('up',0)*n('down',0), #local Hamiltonian
   n_cycles  = 100000, #Monte-Carlo cycles
   length_cycle = 100, #length of a cycle
   n_warmup_cycles = 1000,#thermalization cycles
   n_w_b_vertex = 2, #bosonic freqs. in vertex
   n_w_f_vertex = 20,#fermionic freqs. in vertex
   measure_nn = True, #density-density corr.
   measure_nnw = True, #density-density corr. fct
   measure_gw = True, #one-particle Green fct
   measure_fw = True, #one-particle improved est.
   measure_g2w = True,#3-leg corr.function
   measure_g3w = True,#4-leg corr. function
  )

``h_int`` is the local many-body Hamiltonian, it is an operator expression. ``n_cycles``, ``length_cycle`` and ``n_warmup_cycles`` respectively set the number of configuration measurements, the number of Monte-Carlo updates between two measurements and the number of thermalization cycles. ``n_w_b_vertex`` and ``n_w_f_vertex`` set the number of bosonic and fermionic Matsubara frequencies in the multi-variable correlation functions, and the last six options switch on the measurement of the corresponding observables.

This method comes with mandatory and optional parameters, whose complete list can be found :doc:`here <solver_core>`.

To know more about the solver, visit one of the following pages:

.. toctree::
   :maxdepth: 1

   definitions
   solver_core
   dyn_interactions
   measures
   moves
