.. _measurements:

Measurements 
============

Below is a list of all observables that can be measured by the solver, along with their definitions. Details on the 
implementation of the measurements can be found in the `PhD thesis of T. Ayral <https://hal.archives-ouvertes.fr/tel-01247625>`_ (chapter 11). Each measurement can be 
turned on or off via the corresponding parameter of the ``solve`` method. 

Imaginary time Green's function
*******************************

The imaginary time Green's function is defined as 

.. math::

    G_{ij}^{A}(\tau) = - \langle T_{\tau} c_{Ai}(\tau) c^{\dagger}_{Aj}(0) \rangle, 

where :math:`A` is a block index and :math:`i, j` are inner indices within the block. The block structure 
of :math:`G_{ij}^A(\tau)` (valid values of the indices) is set by ``gf_struct`` in the ``constr_params``.
It is measured on a uniform time grid, whose number of points is set by ``n_tau_G`` in the ``solve_params``
(defaults to ``n_tau`` in ``constr_params``).

.. warning::

    The value of ``n_tau`` supplied in the ``constr_params`` and the number of points in the :math:`\tau` grid of
    the :math:`\Delta(\tau)` input must match. 


The measurement is turned on by setting ``measure_G_tau`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.G_tau`` attribute of the solver object. 

Self-energy improved estimator
******************************

The impurity self-energy can be computed via the Dyson equation 

.. math::

    \Sigma(i\omega_n) = i \omega_n + \mu - G^{-1}(i \omega_n)

from the Green's function :math:`G_{ij}^A(\tau)`. However, this procedure suffers from numerical 
instability at high Matsubara frequencies. A numerically stable way of computing the self-energy is 
provided by measuring the improved estimator 

.. math::

    F_{ij}^A (\tau) = \langle T_{\tau} c_{iA}(\tau) c_{jA}^{\dagger}(0) I_i^{A}(\tau) \rangle.

In the absence of :math:`\mathcal{J}_{\perp}` interactions, 

.. math::

    I_i^A (\tau) = \int_0^{\beta} d\tau' \sum_k \mathcal{U}^A_{ik}(\tau - \tau') n_{Ak}(\tau'). 

The expression of :math:`I_i^A(\tau)` in the presence of :math:`\mathcal{J}_{\perp}` interactions can be found 
in the `PhD thesis of T. Ayral <https://hal.archives-ouvertes.fr/tel-01247625>`_ (Eq. 11.41). 

.. note::

    In the presence of :math:`\mathcal{J}_{\perp}` interactions, the measurement of :math:`F(\tau)` is possible 
    only if the spin-spin interactions are rotationally invariant (:math:`\mathcal{J}_{\perp}(\tau) = \mathcal{J}_z(\tau)`).

Once the improved estimator is known, the self-energy is obtained according to 

.. math::

    \Sigma^A(i\omega_n) = [G^A]^{-1}(i\omega_n) F^A(i \omega_n).

The block structure and imaginary time grid for the improved estimator are the same as for the Green's function. 
The measurement is turned on by setting ``measure_F_tau`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.F_tau`` attribute of the solver object. 

Density
*******

The average density is 

.. math::

    \langle n^A_i \rangle = \frac{1}{\beta} \int_0^{\beta} d \tau \langle n^A_i(\tau) \rangle. 

Here :math:`A` represents a block name and :math:`i` an index within the block. The measurement
is turned on by setting ``measure_densities`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.densities`` attribute of the solver object, as a Python
dictionary that associates a numpy vector to a block name. For example, the average density in the first 
color of the spin up block is accessed as ``results.densities["up"][0]``. 

Static density correlation function
***********************************

The static density correlation function is defined as 

.. math::

    \langle n^A_i n^B_j \rangle = \frac{1}{\beta} \int_0^{\beta} d \tau \langle n^A_i(\tau) n^A_j(\tau) \rangle. 

Here :math:`A, B` represents block names and :math:`i, j` indices within the block. The measurement
is turned on by setting ``measure_nn_static`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.nn_static`` attribute of the solver object, as a Python 
dictionary that associates a matrix to a pair of block names. For example the static density correlation function 
of the spin up block with itself is acccessed as ``results.nn_static[("up", "up")]``. 

Dynamic density correlation function
************************************

The dynamic density correlation function is defined as 

.. math::

    \chi^{AB}_{ij}(\tau) =  \langle T_{\tau} n^A_i(\tau) n^B_j(0) \rangle. 

Here the indices :math:`i, j = 0, \dots N - 1` represents colors (irrespective of the block structure).
:math:`\chi_{ij}(\tau)` is measured on a uniform time grid, whose number of points is set by ``n_tau_chi2`` in the ``solve_params``
(which defaults to ``n_tau_bosonic`` from ``constr_params``). 

The measurement is turned on by setting ``measure_nn_tau`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.nn_tau`` attribute of the solver object, as a 
``Block2Gf``. For example, the correlation function in the first color of the spin up block is accessed as 
``results.nn_tau["up", "up"][0, 0]``. 

Perpendicular spin-spin correlation function
********************************************

The perpendicular spin-spin correlation function is defined as 

.. math::

    \chi^{\perp}(\tau) =  \langle T_{\tau} s^x(\tau) s^x(0) \rangle. 

:math:`\chi^{\perp}(\tau)` is measured on a uniform time grid, whose number of points is set by ``n_tau_chi2`` in the ``solve_params``
(which defaults to ``n_tau_bosonic`` from ``constr_params``).
This measurement is useful if rotational invariance is broken (for instance, in the presence of a Zeeman field). It is 
implemented for a single orbital only. Otherwise, all components of the spin-spin correlation function can be determined from :math:`\chi_{ij}(\tau)`, with better statistics. 

The measurement is turned on by setting ``measure_Sperp_tau`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.Sperp_tau`` attribute of the solver object, as a matrix-valued
``GfImTime`` with size :math:`1 \times 1`.

State histogram
***************

This measurement determines the occupation probabilities of the non-interacting impurity eigenstates. 
Formally, these are the diagonal elements of the impurity density matrix expressed in the occupation
number basis. For example, in the case of an impurity with 2 colors, the eigenstates are 
:math:`|00\rangle, |10\rangle, |01 \rangle, |11\rangle`. 

The measurement is turned on by setting ``measure_state_hist`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.state_hist`` attribute of the solver object, as a numpy array of size
:math:`2^N`. The index of the state :math:`|n_0, n_1, \dots n_N \rangle` in the histogram is given by :math:`\sum_{i = 0}^{N - 1} n_i 2^i`. 

Average sign
************

This measurement computes the average sign of the weight of the configuration. 
It is always on. The result of the accumulation is accessible through the ``results.average_sign`` 
attribute of the solver object as a double precision scalar. 

Perturbation order histograms
*****************************

This measurement determines the histograms of the perturbation orders in :math:`\Delta(\tau)` and :math:`\mathcal{J}_{\perp}(\tau)`. 
The measurement is turned on by setting ``measure_pert_order`` in the ``solve_params`` to ``True``. The results of the 
accumulation are accessible through the ``results.pert_order_Delta`` and ``results.pert_order_Jperp``
attributes of the solver, as TRIQS histogram objects. The average orders can also be directly accessed via 
``results.average_order_Delta`` and ``results.average_order_Jperp``. 
