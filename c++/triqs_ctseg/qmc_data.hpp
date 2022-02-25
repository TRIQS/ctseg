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

    bool const dynamical_U, jperp_interactions, full_spin_rot_inv;

    gf<imtime> const K, Kprime, Kperpprime;

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




    // Construction
    qmc_data(params_t const &p, block_gf_const_view<imtime> delta) delta(map([](gf_const_view<imtime> d) { return real(d); }, delta)) {

      dets.clear();
      for (auto const &bl : range(delta.size())) {
#ifdef HYBRIDISATION_IS_COMPLEX
        dets.emplace_back(delta_block_adaptor(delta[bl]), p.det_init_size);
#else
        if (!is_gf_real(delta[bl], 1e-10)) {
          if (p.verbosity >= 2) {
            std::cerr << "WARNING: The Delta(tau) block number " << bl << " is not real in tau space\n";
            std::cerr << "WARNING: max(Im[Delta(tau)]) = " << max_element(abs(imag(delta[bl].data()))) << "\n";
            std::cerr << "WARNING: Dissregarding the imaginary component in the calculation.\n";
          }
        }
        dets.emplace_back(delta_block_adaptor(real(delta[bl])), p.det_init_size);
#endif
        dets.back().set_singular_threshold(p.det_singular_threshold);
        dets.back().set_n_operations_before_check(p.det_n_operations_before_check);
        dets.back().set_precision_warning(p.det_precision_warning);
        dets.back().set_precision_error(p.det_precision_error);
      }
    }

    K      = gf<imtime>{K_[0].mesh(), make_shape(n_color, n_color)};
    Kprime = gf<imtime>{K_[0].mesh(), make_shape(n_color, n_color)};
    for (auto const &t : K.mesh()) {
      for (int i = 0; i < n_color; i++) {
        auto bl_ind_i = color_to_block_and_inner_index_impl(i, gf_struct);
        for (int j = 0; j < n_color; j++) {
          auto bl_ind_j   = color_to_block_and_inner_index_impl(j, gf_struct);
          int bl_ij       = bl_ind_i.first * gf_struct.size() + bl_ind_j.first;
          K[t](i, j)      = K_[bl_ij](t)(bl_ind_i.second, bl_ind_j.second);
          Kprime[t](i, j) = Kprime_[bl_ij](t)(bl_ind_i.second, bl_ind_j.second);
        } // j
      }   // i
    }     // t
  }
};
} // namespace triqs_ctseg
