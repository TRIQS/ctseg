
Authors
*******

The continuous-time hybridization-expansion segment solver has been written by
T. Ayral, H. Hafermann, P.Delange, M. Ferrero and O. Parcollet.


Useful references
==================

This implementation is based on
a hybridization expansion of the partition function as described in reference
[#werner06]_, with dynamical interactions in the charge and longitudinal spin channel (as described in [#werner07]_) as well as transverse spin channel (as described in [#otsuki13]_).
The computation of the one-particle Green’s
function has been improved with the use of Legendre polynomials following
reference [#boehnke]_. The various improved estimators implemented in the code follow [#hafermann12]_ and [#hafermann14]_.
The code is based on the TRIQS toolbox [#parcollet15]_.

License
=======

The segment solver is published under the `GNU General Public License, version 3
<http://www.gnu.org/licenses/gpl.html>`_.

Quotation
=========

This application is a part of our scientific work and we would appreciate if
projects using it will include a citation to the following relevant papers.  

.. [#werner06] `P. Werner, A. Comanac, L. de' Medici, M. Troyer, and A. J. Millis, Phys. Rev. Lett. 97, 076405 (2006) <http://link.aps.org/doi/10.1103/PhysRevLett.97.076405>`_ 
.. [#werner07] `P. Werner and A. J. Millis, Phys. Rev. B 74, 155107 (2006) <http://link.aps.org/doi/10.1103/PhysRevB.74.155107>`_ 
.. [#otsuki13] `J. Otsuki, Phys. Rev. B 87, 125102 (2013) <http://journals.aps.org/prb/ abstract/10.1103/PhysRevB.87.125102>`_
.. [#boehnke] `L. Boehnke, H. Hafermann, M. Ferrero, F. Lechermann, and O. Parcollet, Phys. Rev. B 84, 075145 (2011) <http://link.aps.org/doi/10.1103/PhysRevB.84.075145>`_ 
.. [#hafermann12] `H. Hafermann, K. R. Patton, P. Werner, Phys. Rev. B 85, 205106 (2012) <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.205106>`_
.. [#hafermann14] `H. Hafermann, Phys. Rev. B 89, 235128 (2014) <http://journals.aps.org/prb/pdf/10.1103/PhysRevB.89.235128>`_
.. [#parcollet15] `O. Parcollet, M. Ferrero, T. Ayral, H. Hafermann, P. Seth, I. S. Krivenko, Comp. Phys. Commun. 196, 398 (2015) <https://doi.org/10.1016/j.cpc.2015.04.023>`_

If you find the application useful, giving proper reference and citation is
indeed a simple way to help convincing funding sources that such projects are
useful for our community and should be supported.

Disclaimer
==========

The program is provided as is, i.e. WITHOUT ANY WARRANTY of any kind, as
stated in the license.  In particular, its authors and contributors will take
no responsability for any possible bugs or any improper use of these programs,
including those resulting in incorrect scientific publications.
