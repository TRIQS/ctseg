/*******************************************************************************
 * CTSEG: TRIQS hybridization-expansion segment solver
 *
 * Copyright (C) 2013-2018 by T. Ayral, H. Hafermann, P. Delange, M. Ferrero, O.
 *Parcollet
 *
 * CTSEG is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * CTSEG is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * CTSEG. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#pragma once
//#include <triqs/utility/first_include.hpp>
#include "./qmc_parameters.hpp"
#include "./types.hpp"

namespace triqs_ctseg {
/// evaluation of the 3 leg vertex
/**

  If one or two of the two-frequency correlation functions have been measured
and the parameter ``evaluate_vertex`` is set to ``true``, the following vertex
function is computed in the charge channel:

.. math::

 \lambda_\sigma(i\omega,i\nu) =
\frac{1}{\chi^{\text{ch}}(i\nu)}\left(\frac{G^{(2),\text{con}}_{\sigma}(i\omega,i\nu)}{G_\sigma(i\omega)G_\sigma(i\omega+i\nu)}
-1\right)

where the charge susceptibility :math:`\chi^\text{ch}(i\nu)` is defined as

.. math::
   \chi^\text{ch}(i\nu) = -\sum_{\sigmab} \left[\chi_{\sigmab}(i\nu)-\langle
n_\sigma\rangle \langle n_b\rangle \delta_{\nu} \right]

and where, depending on which two-frequency correlation functions have been
measured, the connected part is computed in either of the following ways:

.. math::

 \begin{align}
 G^{(2),\text{con}}_{\sigma}(i\omega,i\nu) &= \sum_b
G^{(2)}_{\sigmab}(i\omega,i\nu) - \sum_b
G^{(2),\text{disc}}_{ab}(i\omega,i\nu)\\
 G^{(2),\text{con}}_{\sigma}(i\omega,i\nu) &= G_\sigma(i\omega) \sum_b
F_{\sigmab}^{(2)}(i\omega,i\nu)-F_\sigma(i\omega) \sum_b
G_{\sigmab}^{(2)}(i\omega,i\nu)\\ G^{(2),\text{con}}_{\sigma}(i\omega,i\nu) &=
\Big[G_\sigma(i\omega) \sum_b F_{\sigmab}^{(2)}(i\omega,i\nu)-F_\sigma(i\omega)
\sum_b G_{\sigmab}^{(2),\text{disc}}(i\omega,i\nu)\Big]/[1+F_\sigma(i\omega)].
 \end{align}

Here the disconnected part of the correlation function has been defined as

.. math::

 G^{(2),\text{disc}}_{\sigmab}(i\omega,i\nu) = \beta G_\sigma(i\omega)\langle
n_b\rangle\delta_{\nu}-G_\sigma(i\omega) G_\sigma(i\omega+i\nu)\delta_{\sigmab}.

  */
/*
std::pair<array<gf_2w_t, 1> ,array<gf_2w_t, 1> > evaluate_2w_vertex(const
mpi::communicator &c, const qmc_parameters *params_, gf<block_index, gf<imfreq>
> &gw, gf<block_index, gf<imfreq> > &fw, block_f_om_nu_tv_t &g2w,
block_f_om_nu_tv_t &f2w, gf<imfreq> &nnw, matrix<double> &nn_matrix);
                       */
/// Evaluation of the 4-leg vertex for the 4-point correlation function
/**
 *If one or two of the three-frequency correlation functions have been measured
 and the parameter ``evaluate_vertex`` is set to ``true``, the following vertex
 function is computed at the end of the simulation:

  $$\gamma_{\sigma\sigma'}(i\omega,i\omega',i\nu) =
 \frac{G^{2,\text{con}}_{\sigma\sigma'}(i\omega,i\omega',i\nu)}{G_\sigma(i\omega)G_\sigma(i\omega+i\nu)G_\sigma'(i\omega'+i\nu)G_\sigma'(i\omega')}$$

 *Depending on which two-frequency correlation functions have been measured, the
 connected part is computed in either of the following ways:


  $$G^{2,\text{con}}_{\sigma}(i\omega,i\omega',i\nu) =
 G^{2}_{\sigma\sigma'}(i\omega,i\omega',i\nu) -
 G^{2,\text{disc}}_{\sigma\sigma'}(i\omega,i\omega',i\nu)$$

  $$G^{2,\text{con}}_{\sigma}(i\omega,i\omega',i\nu) = G_\sigma(i\omega)
 F_{\sigma\sigma'}^{2}(i\omega,i\omega',i\nu)-F_\sigma(i\omega)
 G_{\sigma\sigma'}^{2}(i\omega,i\omega',i\nu)$$

  $$G^{2,\text{con}}_{\sigma}(i\omega,i\omega',i\nu) = \Big[G_\sigma(i\omega)
 F_{\sigma\sigma'}^{2}(i\omega,i\omega',i\nu)-F_\sigma(i\omega)
 G_{\sigma}^{2,\text{disc}}(i\omega,i\omega',i\nu)\Big]/[1+F_\sigma(i\omega)].$$

 *The disconnected part of the correlation function has been defined as

  $$G^{2,\text{disc}}_{\sigma\sigma'}(i\omega,i\omega',i\nu) = \beta
 G_{\sigma}(i\omega)G_\sigma'(i\omega')\delta_{\nu}-\beta G_\sigma(i\omega)
 G_\sigma(i\omega+i\nu)\delta_{\omega,\omega'}\delta_{\sigma\sigma'}.$$

  */
void evaluate_3w_vertex(block_gf<imfreq> const &gw, block_gf<imfreq> const &fw,
                        gf_3w_container_t const &g3w,
                        gf_3w_container_t const &f3w, bool measure_g3w,
                        bool measure_f3w, std::string fname);
} // namespace triqs_ctseg
