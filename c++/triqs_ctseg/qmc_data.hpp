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
#include "types.hpp"
#include "params.hpp"
#include "inputs_t.hpp"

namespace triqs_ctseg {

  using time_point_factory_t = triqs::utility::time_segment;

  /// A lambda to adapt the Delta function for the call of the det.
  struct delta_block_adaptor {
    gf<imtime, matrix_real_valued> delta; // make a copy. Needed in the real case anyway.

    double operator()(std::pair<qmc_time_pt, int> const &x, std::pair<qmc_time_pt, int> const &y) const {
      //det_scalar_t res = delta[closest_mesh_pt(double(x.first - y.first))](x.second, y.second);
      double res = delta(double(x.first - y.first))(x.second, y.second);
      return (x.first >= y.first ? res : -res); // x,y first are time_pt, wrapping is automatic in the - operation, but need to
                                                // compute the sign
    }
  };

  // ---------------------------------------------------
  /// Working data
  struct work_data_t {
    work_data_t(param_t const &params, inputs_t const &inputs);

    int n_color;
    double beta;
    nda::vector<double> mu;
    nda::matrix<double> U;

    bool has_Dt, has_jperp;
    gf<imtime> K, Kprime;

    // FIXME off diagonal delta ??
    block_gf<imtime, delta_target_t> delta; // Hybridization function
    using det_t = det_manip::det_manip<delta_block_adaptor>;
    std::vector<det_t> dets; // The determinants
  };

} // namespace triqs_ctseg
