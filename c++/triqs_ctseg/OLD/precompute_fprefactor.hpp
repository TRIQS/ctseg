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
#include "./configuration.hpp"
#include "./qmc_parameters.hpp"

namespace triqs_ctseg {

/// Object which precomputes the prefactor $I_{o_1}(\tau)$ for the improved
/// estimator
/**
 *$I_{o_1}$ is used in the computation of the improved estimator
 *$F_{o_1}(\tau)$. It is defined as:
 *  $$I_{o_1}(\tau) = \int_0^\beta d\bar{\tau} \sum_{o_2}
 *\hat{\mathcal{U}}_{o_1,o_2}(\tau-\bar{\tau}) n_{o_2}(\bar{\tau})$$ This object
 *models the concept of mc_generic auxiliary precomputation.
 * @note in doc/notes for more details
 * @note In the measurement all times are accessed for a given orbital and not
 *vice versa; Storage in vector of maps instead of map of vectors keeps the
 *individual maps small
 */
struct precompute_fprefactor {

  const qmc_parameters *params;
  const configuration *config;

  using fprefactor_t = double;
  using fprefactormap_t =
      std::map<qmc_time_t, fprefactor_t, std::greater<qmc_time_t>>;

  std::vector<fprefactormap_t> fprefactors;

  precompute_fprefactor(const qmc_parameters *params_,
                        const configuration *config_);

  /// accessor
  fprefactor_t get(int color, qmc_time_t t) {
    try {
      return fprefactors[color].at(t);
    } catch (const std::out_of_range &exc) {
      TRIQS_RUNTIME_ERROR << "operator not found in fprefactors: " << exc.what()
                          << "\n config = " << *config << "\n t = " << t
                          << "\t c = " << color;
    }
  }

  /// call operator which performs the actual computation
  void operator()();
};

} // namespace triqs_ctseg
