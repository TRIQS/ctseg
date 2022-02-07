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
#include "../params.hpp"
#include "../precompute_fprefactor.hpp"
#include "../qmc_parameters.hpp"

namespace triqs_ctseg {

using namespace triqs::gfs;
using namespace triqs::mesh;

/// measures Fourier transform of $G(\tau)$ and $F(\tau)$
/**
 *Measurement of
  $$X^\sigma_a(i\omega_n) = \int_0^\beta d\tau\, e^{i\omega_n\tau}
 X^\sigma_a(\tau)$$
 *
 * with $X=G,F$ defined in [[measure_gt]].
 *
 *The switches for these measurements are ``measure_gw`` and ``measure_fw``. The
 number of frequencies is specified through ``n_w``.
 *
 *This implementation uses [[precompute_fprefactor]].
 *
 *The self-energy is automatically computed at the end of the simulations as
 $$\Sigma_{ab}^{\sigma}(i\omega_{n})=\sum_c
 F_{ac}^{\sigma}(i\omega_{n})\left[G^{-1}\right]_{cb}^{\sigma}(i\omega_{n})$$
 *
 * if the measurement of the improved estimator is turned on, or
 $$\Sigma_{ab}^{\sigma}(i\omega_{n})\equiv\left[\mathcal{G}^{-1}\right]_{ab}^{\sigma}(i\omega_{n})-\left[G^{-1}\right]_{ab}^{\sigma}(i\omega_{n})$$
 *
 * otherwise.
 */
struct measure_gw {

  const qmc_parameters *params;
  const solve_params_t *p;
  const configuration *config;

  block_gf<imfreq> &gw;
  block_gf<imfreq> &fw;
  block_gf<imfreq> &sigmaw;
  block_gf_const_view<imfreq> g0w;

  std::shared_ptr<precompute_fprefactor> fprefactor;

  int n_w;
  double beta;
  double Z;

  vector<dcomplex> c_exp, cdag_exp;
  vector<int> c_inner_index, cdag_inner_index;

  /// constructor
  measure_gw(const qmc_parameters *params_, const solve_params_t *p_,
             const configuration *config_,
             std::shared_ptr<precompute_fprefactor> fprefactor_,
             block_gf<imfreq> &gw_, block_gf<imfreq> &fw_,
             block_gf<imfreq> &sigmaw_, block_gf_const_view<imfreq> g0w_);

  /// accumulate the Green's function
  void accumulate(double s);

  /// reduce and normalize G
  void collect_results(mpi::communicator const &c);
};
} // namespace triqs_ctseg
