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
#include "../configuration.hpp"
#include "../precompute_Mw.hpp"
#include "../precompute_fprefactor.hpp"
#include "../precompute_nw.hpp"
#include "../qmc_parameters.hpp"
#include "../types.hpp"

namespace triqs_ctseg {

/// Measures Fourier transform of  $\chi^{\sigma,\sigma'}_{ab}(\tau,\tau')$ and
/// $F^{2,\sigma,\sigma'}_{ab}(\tau,\tau')$
/**
 *This class implements the measurement of:

 $$X_{ab}(i\omega_n,i\nu_m) = \int_0^\beta d\tau \int_0^\beta d\tau'
e^{i\omega_n\tau} e^{i\nu_m\tau'} X_{ab}(\tau,\tau')$$

 with $X=G^{2},F^{2}$ defined as:

 $$G^{2,\sigma\sigma'}_{abc}(\tau,\tau') = -\langle T_\tau
c_{a\sigma}(\tau)c_{b\sigma}^\dagger(0)n_{c\sigma'}(\tau') \rangle$$

 * and

 $$F^{2,\sigma\sigma'}_{abc}(\tau,\tau') = -\int_0^\beta d\tilde{\tau}
\sum_{d\bar{\sigma}} \langle T_\tau n_{d\bar{\sigma}}(\tilde{\tau})
\mathcal{U}^{\sigma\bar{\sigma}}_{ad}(\tilde{\tau}-\tau)
c_{a\sigma}(\tau)c_{b\sigma}^\dagger(0)n_{c\sigma'}(\tau') \rangle$$

* The number of fermionic (bosonic) frequencies is specified through the
parameters ``n_w_f_vertex`` (``n_w_b_vertex``).
* This implementation uses [[precompute_Mw]].
 @note Currently only the frequency measurement is provided.
*/
struct measure_g2w {

  const qmc_parameters *params;
  const configuration *config;

  block_f_om_nu_tv_t &g2w;
  block_f_om_nu_tv_t &f2w;

  int n_w_fermionic, n_w_bosonic;
  double beta, Z;

  std::shared_ptr<precompute_Mw> Mw;
  std::shared_ptr<precompute_nw> nw;

  accumulator<double> g2w_0_0_stack   = {0.0, -1, -1};
  accumulator<double> g2w_10_0_stack  = {0.0, -1, -1};
  accumulator<double> g2w_m10_0_stack = {0.0, -1, -1};

  measure_g2w(const qmc_parameters *params_, const configuration *config_,
              std::shared_ptr<precompute_Mw> Mw_,
              std::shared_ptr<precompute_nw> nw_, block_f_om_nu_tv_t &g2w_,
              block_f_om_nu_tv_t &f2w_);

  /// accumulate the Green's function
  /**
    this measures the Fourier transform of $-<T
    c_1(\tau_1)c^*_1(\tau_2)n_2(\tau_3)>$ and  $- \sum_1 <T n_1'(\tau_1) U_1'1
    c_1(\tau_1)c^*_1(\tau_2)n_2(\tau_3)>$
   */
  void accumulate(double s);

  /// reduce and normalize G
  void collect_results(mpi::communicator const &c);
};
} // namespace triqs_ctseg
