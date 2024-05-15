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
It is measured on a uniform time grid, whose number of points is set by ``n_tau`` in the ``constr_params``. 

.. warning::

    The value of ``n_tau`` supplied in the ``constr_params`` and the number of points in the :math:`\tau` grid of
    the :math:`\Delta(\tau)` input must match. 


The measurement is turned on by setting ``measure_gt`` in the ``solve_params`` to ``True``. The result of the 
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

    F_{ij}^A (\tau) = \langle T_{\tau} c_{iA}(\tau) c_{jA}^{\dagger}(0) I_i^{A}(\tau). 

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
The measurement is turned on by setting ``measure_ft`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.F_tau`` attribute of the solver object. 

Density
*******

The average density is 

.. math::

    \langle n_i \rangle = \frac{1}{\beta} \int_0^{\beta} d \tau \langle n_i(\tau) \rangle. 

Here the index :math:`i = 0, \dots N - 1` represents the color (irrespective of the block structure). The measurement
is turned on by setting ``measure_n`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.densities`` attribute of the solver object, as a numpy 
array of size :math:`N`. 

Static density correlation function
***********************************

The static density correlation function is defined as 

.. math::

    \langle n_i n_j \rangle = \frac{1}{\beta} \int_0^{\beta} d \tau \langle n_i(\tau) n_j(\tau) \rangle. 

Here the indices :math:`i, j = 0, \dots N - 1` represents colors (irrespective of the block structure). The measurement
is turned on by setting ``measure_nn`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.nn_static`` attribute of the solver object, as a numpy 
array of size :math:`N \times N`. 

Dynamic density correlation function
************************************

The dynamic density correlation function is defined as 

.. math::

    \chi_{ij}(\tau) =  \langle T_{\tau} n_i(\tau) n_j(0) \rangle. 

Here the indices :math:`i, j = 0, \dots N - 1` represents colors (irrespective of the block structure).
:math:`\chi_{ij}(\tau)` is measured on a uniform time grid, whose number of points is set by ``n_tau_k`` in the ``constr_params``. 

.. warning::

    The value of ``n_tau_k`` supplied in the ``constr_params`` and the number of points in the :math:`\tau` grids of
    the :math:`D(\tau)` and :math:`J_{\perp}(\tau)` inputs must match. 

The measurement is turned on by setting ``measure_nnt`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.nn_tau`` attribute of the solver object, as a matrix-valued
``GfImTime`` with size :math:`N \times N`. 

Perpendicular spin-spin correlation function
********************************************

The perpendicular spin-spin correlation function is defined as 

.. math::

    \chi^{\perp}(\tau) =  \langle T_{\tau} s^x(\tau) s^x(0) \rangle. 

:math:`\chi^{\perp}(\tau)` is measured on a uniform time grid, whose number of points is set by ``n_tau_k`` in the ``constr_params``. 
This measurement is useful if rotational invariance is broken (for instance, in the presence of a Zeeman field). Otherwise, 
all components of the spin-spin correlation function can be determined from :math:`\chi_{ij}(\tau)`, with better statistics. 

The measurement is turned on by setting ``measure_sperpt`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.sperp_tau`` attribute of the solver object, as a matrix-valued
``GfImTime`` with size :math:`1 \times 1`.

State histogram
***************

This measurement determines the occupation probabilities of the non-interacting impurity eigenstates. 
Formally, these are the diagonal elements of the impurity density matrix expressed in the occupation
number basis. For example, in the case of an impurity with 2 colors, the eigenstates are 
:math:`|00\rangle, |10\rangle, |01 \rangle, |11\rangle`. 

The measurement is turned on by setting ``measure_statehist`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.state_hist`` attribute of the solver object, as a numpy array of size
:math:`2^N`. The index of the state :math:`|n_0, n_1, \dots n_N \rangle` in the histogram is given by :math:`\sum_{i = 0}^{N - 1} n_i 2^i`. 

Average sign
************

This measurement computes the average sign of the weight of the configuration. 
The measurement is turned on by setting ``measure_sign`` in the ``solve_params`` to ``True``. The result of the 
accumulation is accessible through the ``results.sign`` attribute of the solver object as a double precision scalar. 

Perturbation order histograms
*****************************

This measurement determines the histograms of the perturbation orders in :math:`\Delta(\tau)` and :math:`\mathcal{J}_{\perp}(\tau)`. 
The measurement is turned on by setting ``measure_statehist`` in the ``solve_params`` to ``True``. The results of the 
accumulation are accessible through the ``results.perturbation_order_histo_Delta`` and ``results.perturbation_order_histo_Jperp``
attributes of the solver, as TRIQS histogram objects. 
