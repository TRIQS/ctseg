.. module:: ctseg

.. _aim:

The Anderson Impurity Model with dynamical interactions
========================================================

Definition
-----------

The Anderson impurity model (AIM) is defined by the action:



.. math::

  S &= \iint_0^\beta \mathrm{d}\tau \mathrm{d}\tau' \sum_{ab\sigma} c^{*}_{a\sigma}(\tau) \left( - \left[\mathcal{G}^{-1}_0\right]_{ab}^\sigma(\tau-\tau') \right) c_{b\sigma}(\tau')\\ &\;\;+ \frac{1}{2}\iint_{0}^{\beta}\mathrm{d}\tau\mathrm{d}\tau'\sum_{ab}\sum_{\sigma\sigma'}\mathcal{U}_{ab}^{\sigma\sigma'}(\tau-\tau')n_{a\sigma}(\tau)n_{b\sigma'}(\tau') \\ &\;\;+ \frac{1}{2}\iint_{0}^{\beta}\mathrm{d}\tau\mathrm{d}\tau'\sum_{a}\sum_{\xi=x,y,z}s_{a}^{\xi}(\tau)\mathcal{J}_{a}^{\xi}(\tau-\tau')s_{a}^{\xi}(\tau')

where :math:`\mathcal{G}_0`  is the non-interacting Green's function, related to the hybridization function :math:`\Delta` by the relation :

.. math::

  [\mathcal{G}^{-1}_{0}]_{ab}^\sigma (i \omega_n) = (i \omega_n + \mu_{a\sigma})\delta_{ab} - \Delta_{ab}^\sigma(i \omega_n).


:math:`a, b \dots` are the orbital indices, :math:`n_{a\sigma}\equiv c^\dagger_{a\sigma}c_{a\sigma}`, :math:`s_{a}^{\xi}\equiv\frac{1}{2}\sum_{\sigma\sigma'}c_{a\sigma}^{*}\sigma_{\sigma\sigma'}^{\xi}c_{a\sigma'}`, :math:`\mathcal{U}(\tau)` is the density-density dynamical interaction and :math:`\mathcal{J}_\perp` is the transverse part of the dynamical spin interactions. :math:`\mathcal{U}` is in general composed of a static part :math:`U` and of a retarded part :math:`\mathcal{D}_0(\tau)`:

.. math::

  \mathcal{U}(\tau) = U\delta(\tau) + \mathcal{D}_0(\tau)


Hybridization expansion and segment picture
--------------------------------------------
Solving the Anderson impurity model amounts to computing its correlation functions like the one-particle Green's function

.. math::

  G^\sigma_{ab}(\tau) \equiv - \langle T c_{a\sigma}(\tau) c^{\dagger}_{b\sigma}(0) \rangle

In order to do so, one can expand one of the terms of the Anderson action and numerically evaluate the corresponding Feynman diagrams. The hybridization expansion consists in expanding the partition function in powers of the hybridization function :math:`\Delta(\tau)`. In the case of density-density interactions :math:`U_{ab} n_{a}(\tau) n_{b}(\tau)` the diagrams one obtains can be represented by sets of occupied or empty segments, hence the name "segment picture". This segment picture allows for an analytical formula for and hence fast computation of the trace of the impurity operators :


.. math::

  w_U = \exp\left(-\sum_{a,b,b<a} U_{ab} l^{\mathrm{overlap}}_{ab} + \mu \sum_a l_a\right)



Density-density dynamical interactions
----------------------------------------
Dynamical interactions do not alter the segment picture. The trace can still be computed without great numerical expense. Dynamical interactions just result in a shift of the static :math:`U_{ab}` by :math:`-2 K'_{ab}(0)` and of the chemical potential :math:`\mu_a`  of :math:`K'_a(0)` , as well as an additional term in the trace:

.. math::

  w_{\mathcal{D}_0} = \exp\left(\sum_{ab} \sum_{\mathrm{op} i_a,j_b} s_{i_{a}} s_{j_{b}} K_{ab}(\tau_{i_a} - \tau_{j_b}) \right)


where :math:`s_{i_a} = 1`  if the operator at time :math:`\tau_{i_a}` is a creation operator,  :math:`s_{i_a} = -1` otherwise.

The kernel :math:`K(\tau)` is defined as 

.. math::

  K''(\tau) = \mathcal{D}_0(\tau) \\
  K(0) = K(\beta) = 0

It is obtained from the retarded interaction :math:`\mathcal{D}_0(\tau)` through the expression:


.. math::

  K(\tau) = \frac{2}{\beta} \sum_{n>0} \frac{\mathcal{D}_0(i\nu_n)-\mathcal{D}_0(0)}{\nu_n^2} \left(  1-\cos(\nu_n \tau) \right)

:math:`\mathcal{D}_0(i\nu_n)` is set in the solver using the accessor `D0_iw`.

Transverse spin-spin dynamical interactions
---------------------------------------------

Contrary to density-density dynamical interactions, transverse spin-spin dynamical interaction (terms :math:`\xi=x,y` of the AIM action above) necessitate an additional expansion in powers of :math:`\mathcal{J}_\perp`. 

:math:`\mathcal{J}_\perp(i\nu_n)` is set in the solver using the accessor `Jperp_iw`.

.. note::

   In this implementation, transverse spin-spin dynamical interactions are limited to the one-band case
