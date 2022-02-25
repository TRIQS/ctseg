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
#include "./params.hpp"
//#include "./util.hpp"
#include <triqs/utility/time_pt.hpp>

namespace triqs_ctseg {

  using time_point_factory_t = triqs::utility::time_segment;

  /// Parameters for QMC
  struct qmc_data_t {

    int const n_color;
    double const beta;
    nda::matrix<double> const U;
    nda::vector<double> const mu;

    bool const dynamical_U, jperp_interactions, full_spin_rot_inv; // FIXME: optional Jperp? 

    gf<imtime> const K, Kprime;

    // FIXME off diagonal delta ??
    block_gf<imtime, delta_target_t> delta; // Hybridization function

    /// A lambda to adapt the Delta function for the call of the det.
    struct delta_block_adaptor {
      gf<imtime, delta_target_t> delta_block; // make a copy. Needed in the real case anyway.

      det_scalar_t operator()(std::pair<qmc_time_pt, int> const &x, std::pair<qmc_time_pt, int> const &y) const {
        det_scalar_t res = delta_block[closest_mesh_pt(double(x.first - y.first))](x.second, y.second);
        return (x.first >= y.first ? res : -res); // x,y first are time_pt, wrapping is automatic in the - operation, but need to
                                                  // compute the sign
      }
    };

    std::vector<det_manip::det_manip<delta_block_adaptor>> dets; // The determinants

};
} // namespace triqs_ctseg
