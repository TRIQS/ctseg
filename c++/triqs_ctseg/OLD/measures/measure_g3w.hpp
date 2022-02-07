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
#include "../qmc_parameters.hpp"
#include "../types.hpp"

namespace triqs_ctseg {

/// This measures $\chi_{abcd}^{\sigma\sigma'}(i\omega_1,i\omega_2,i\omega_3)$
/// and $F^{3,\sigma\sigma'}_{abcd}(i\omega_1,i\omega_2,i\omega_3)$
/**
 *The four-point correlation function is defined as:
 *
 *$$\chi^{\sigma\sigma'}_{abcd}(i\omega, i\omega',i\Omega) = G^{2,\sigma,
 *\sigma'}_{abcd}(i\omega, i\omega',i\Omega) = \langle c_{a\sigma}(i\omega)
 *c^\dagger_{b\sigma}(i\omega+i\Omega) c_{c\sigma'}(i\omega'+i\Omega)
 *c^\dagger_{d\sigma'}(i\omega') \rangle$$
 *
 * Its improved estimator is the Fourier transform of
 *
 *$$F^{3,\sigma\sigma'}_{abcd}(\tau,\tau',\tau'') = \int \mathrm{d}\bar{\tau}
 *\sum_{e\bar{\sigma}} \mathcal{U}^{\sigma\bar{\sigma}}_{ae}(\bar{\tau}-\tau)
 *\langle c_{a\sigma}(\tau) c^\dagger_{b\sigma}(\tau') c_{c\sigma'}(\tau'')
 *c^\dagger_{d\sigma'}(0) \rangle$$
 *
 * The vertex corresponding to this correlation function is evaluated separately
 *(see [[evaluate_3w_vertex]]).
 * @warning The result is reduced only on the master node.
 */
struct measure_g3w {

  const qmc_parameters *params;
  const configuration *config;

  gf_3w_container_t &g3w;
  gf_3w_container_t &f3w;

  int n_w_fermionic, n_w_bosonic;
  double beta, Z;

  std::shared_ptr<precompute_Mw> Mw;

  /// constructor
  measure_g3w(const qmc_parameters *params_, const configuration *config_,
              std::shared_ptr<precompute_Mw> Mw_, gf_3w_container_t &g3w,
              gf_3w_container_t &f3w_);

  /// accumulation
  void accumulate(double s);

  /// reduce and normalize G
  void collect_results(mpi::communicator const &c);
};

} // namespace triqs_ctseg
