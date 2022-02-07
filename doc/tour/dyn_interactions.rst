.. _documentation:

Example with dynamical charge and spin interactions
=========================================================

Below is a short example demonstrating how to set up a computation in the presence of dynamical interactions in the charge and longitudinal and transverse spin channels.


.. runblock:: python

  from triqs.gf import *
  import triqs.utility.mpi as mpi
  from triqs_ctseg import Solver

  gf_struct = [("up",[0]), ("down",[0])]
  S = Solver(beta = 20.0,  gf_struct = gf_struct)
  #set noninteracting Green's function
  mu, eps0, V = 2.0, 0.3, 1.0
  S.G0_iw << inverse(iOmega_n + mu - V**2 * inverse(iOmega_n - eps0))

  #set dynamical charge and longitudinal spin interactions
  w_ch, w_sp = 1.0, 0.8 #screening frequencies
  g_ch, g_sp = 0.4, 0.7 #coupling strengths
  sigz = lambda s1, s2 : 1 if s1==s2 else -1
  import itertools
  from triqs.gf.descriptors import Function
  for b1, b2 in itertools.product(dict(gf_struct).keys(), repeat=2):
   S.D0_iw[b1+"|"+b2] << Function(lambda w : g_ch**2*2*w_ch/(w**2-w_ch**2)\
                                           + g_sp**2*2*w_sp/(w**2-w_sp**2)*sigz(b1,b2) )
  #set transverse spin interaction
  S.Jperp_iw << Function(lambda w : 4 * g_sp**2 * 2*w_sp/(w**2-w_sp**2))

  from triqs.operators import *
  U = 2.0
  S.solve(h_int = U * n('up',0)*n('down',0),
          n_cycles = 10000, length_cycle = 100, n_warmup_cycles = 100)

  from h5 import *
  if mpi.is_master_node():
   with HDFArchive("dyn_interactions.h5",'w') as A: 
    A['G_tau'] = S.G_tau



Let us decompose the various steps:

*Import the* ``Solver`` *class*::

  from triqs_ctseg import Solver 
 
*Declare and construct a* ``Solver`` *instance*::

  gf_struct = [("up",[0]), ("down",[0])]
  S = Solver(beta = 20.0,  gf_struct = gf_struct)

Here, ``gf_struct`` is a Python dictionary that specifies the block structure of the Green's function, here one :math:`1\times 1` (``[0]`` is of length 1) for each spin (``up``, ``down``).

Some construction parameters are mandatory, others are optional. See :doc:`here <solver_core>` for a complete list.
 
*Initialize the inputs of the solver*. In this script, we want to define the noninteracting Green's function as:

.. math::

  \begin{equation}
  \left[\mathcal{G}_{0}\right]^{\sigma}(i\omega_{n})=\left[i\omega_{n}-\varepsilon_{a\sigma}-\Delta_{ab}^{\sigma}(i\omega_{n})\right]^{-1},\label{eq:G0_def}
  \end{equation}

with 

.. math::

  \Delta^\sigma(i\omega_{n}) =\frac{V^{2}}{i\omega_{n}-\varepsilon_{0}}\\

In the script, this is done via an *accessor* called ``G0_iw``. It is set by the user in the following way::

  S.G0_iw << inverse(iOmega_n + mu - V**2 * inverse(iOmega_n - eps0))

Similarly, we want to define dynamical interactions:

.. math::

  \begin{align}
  \mathcal{D}_{0}^{\sigma\sigma'}(i\Omega) & =\frac{2g_{\mathrm{ch}}^{2}\omega_{\mathrm{ch}}}{(i\Omega)^{2}-\omega_{\mathrm{ch}}^{2}}+\left(-\right)^{\sigma\sigma'}\frac{2g_{\mathrm{sp}}^{2}\omega_{\mathrm{sp}}}{(i\Omega)^{2}-\omega_{\mathrm{sp}}^{2}},\label{eq:D0_example}\\
  \mathcal{J}_{\perp}(i\Omega) & =4\cdot\frac{2g_{\mathrm{sp}}^{2}\omega_{\mathrm{sp}}}{(i\Omega)^{2}-\omega_{\mathrm{sp}}^{2}}.\label{eq:Jperp}
  \end{align}


This is done through the ``D0_iw`` and ``Jperp_iw`` accessors::

  for b1, b2 in itertools.product(dict(gf_struct).keys(), repeat=2):
   S.D0_iw[b1+"|"+b2] << Function(lambda w : g_ch**2*2*w_ch/(w**2-w_ch**2)\
                                           + g_sp**2*2*w_sp/(w**2-w_sp**2)*sigz(b1,b2) )

  S.Jperp_iw << Function(lambda w : 4 * g_sp**2 * 2*w_sp/(w**2-w_sp**2))
 
For a complete list of the accessors (both input and output), see :doc:`here <solver_core>`.

*Call the* ``solve()`` *method*::

  S.solve(h_int = U * n('up',0)*n('down',0),
          n_cycles = 10000, length_cycle = 100, n_warmup_cycles = 100)

``h_int`` is the local many-body Hamiltonian, it is an operator expression. ``n_cycles``, ``length_cycle`` and ``n_warmup_cycles`` respectively set the number of configuration measurements, the number of Monte-Carlo updates between two measurements and the number of thermalization cycles.

The fact that ``Jperp_iw`` is nonzero is automatically detected by the code and the corresponding Monte-Carlo updates are turned on. 

The ``solve()`` method comes with mandatory and optional parameters, whose complete list can be found :doc:`here <solver_core>`.

To know more about the solver, visit one of the following pages:

.. toctree::
   :maxdepth: 1

   definitions
   solver_core
   measures
