.. _welcome:

The Segment Picture Solver
**************************

The :ref:`TRIQS-based <triqslibs:welcome>` hybridization-expansion **segment picture** solver (CTSEG)
can tackle the generic problem of a quantum impurity coupled to an external environement (bath). The "impurity" can
be any set of orbitals, on one or several atoms. The CTSEG solver supports (possibly retarded) density-density
and spin-spin interactions on the impurity. Under these restrictions, it provides better performance than the generic
`CTHYB <https://triqs.github.io/cthyb/latest/>`_ solver, that supports generic local interaction vertices. 
The imaginary time action solved by CTSEG is of the form

.. math::

  \begin{split}
  \mathcal{S}  &= \iint_0^{\beta} \mathrm{d} \tau \mathrm{d} \tau' \sum_{a,b} \left\{ \overline{c}_{a\sigma} (\tau)
  \left( (\partial_{\tau} + \epsilon_{a\sigma})\delta_{ab}^{\sigma \sigma'} \delta_{\tau - \tau'} + \Delta_{ab}^{\sigma \sigma'}(\tau - \tau')\right)
  c_{b\sigma'}(\tau') \right\} \\
  &+ \frac{1}{2} \iint_0^{\beta} \mathrm{d} \tau \mathrm{d} \tau' \sum_{a,b} \mathcal{U}_{ab}(\tau - \tau') n_a(\tau) n_b(\tau') 
  + \frac{1}{2} \iint_0^{\beta} \mathrm{d} \tau \mathrm{d} \tau' \sum_{a, \xi = x, y, z} s_a^{\xi}(\tau) \mathcal{J}_a^{\xi}(\tau - \tau') s_a^{\xi} (\tau')
  \end{split}

.. sidebar:: CTSEG |PROJECT_VERSION|

   This is the homepage of CTSEG |PROJECT_VERSION|.
   For changes see the :ref:`changelog page <changelog>`.
      
      .. image:: _static/logo_github.png
         :width: 65%
         :align: center
         :target: https://github.com/triqs/ctseg

Here :math:`\beta` is the inverse temperature, :math:`a` denote orbital indices, :math:`\sigma` spin indices (:math:`\sigma = \uparrow, \downarrow`),
:math:`n_a \equiv \sum_{\sigma} n_{a\sigma}`, :math:`s_a^{\xi} \equiv \frac{1}{2} \sum_{\sigma \sigma'} \overline{c}_{a\sigma}
\sigma_{\sigma \sigma'}^{\xi} c_{a \sigma'}` and :math:`\sigma^{\xi}` are the Pauli matrices. :math:`\overline{c}_{a\sigma}(\tau)`
and :math:`c_{a\sigma}(\tau)` are the :math:`\beta`-antiperiodic Grassman fields corresponding to the fermion
creation and annihilation operators on the impurity, respectively. :math:`\Delta_{ab}^{\sigma \sigma'}(\tau)` 
is the hybridization function, that accounts for particle exchange between the impurity and the bath, and 
:math:`\mathcal{U}_{ab} (\tau)` and  :math:`\mathcal{J}_{a}^{\xi} (\tau)` are the (dynamical)
density-density and spin-spin interactions, respectively. 

The CTSEG solver carries out a double expansion in the hybridization term and in the perpendicular spin-spin
interaction term to obtain the fully interacting impurity Green's function :math:`G(\tau)` and a range of
other observables. Learn how to use it in the :ref:`documentation`.

.. image:: _static/logo_flatiron.png
   :width: 50%
   :target: "https://www.simonsfoundation.org/flatiron/center-for-computational-quantum-physics/

    
.. toctree::
   :maxdepth: 2
   :hidden:

   install
   documentation
   issues
   ChangeLog.md
   about
