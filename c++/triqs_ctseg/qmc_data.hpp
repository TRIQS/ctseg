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

 // CLEAN THIS ....



using triqs::utility::time_segment;
/// Parameters for QMC
struct qmc_data {

  //
  solve_params_t const &params;


  // data.
  double beta;
  time_segment tau_seg;
  matrix<double> U;
  vector<double> mu;
  bool dynamical_U, jperp_interactions, full_spin_rot_inv;
  gf<imtime> K, Kprime, Kperpprime;
  gf_struct_t gf_struct;
  
  
  int n_color; // OK

  int n_w_f_vertex, n_w_b_vertex;
  
  std::string fname_gammaw;

  qmc_data(int n_color_, double beta_, array<double, 2> U_,
                 array<double, 1> mu_, block_gf<imtime> K_,
                 block_gf<imtime> Kprime_, gf<imtime> Kperpprime_,
                 gf_struct_t const &gf_struct_, bool dynamical_U_,
                 bool jperp_interactions_, bool full_spin_rot_inv_,
                 solve_params_t const &p)
      : beta(beta_), tau_seg(beta), U(U_), mu(mu_), dynamical_U(dynamical_U_),
        jperp_interactions(jperp_interactions_),
        full_spin_rot_inv(full_spin_rot_inv_), Kperpprime(Kperpprime_),
        gf_struct(gf_struct_), n_color(n_color_),
        measure_g2w(p.measure_g2w), measure_f2w(p.measure_f2w),
        measure_g3w(p.measure_g3w), measure_f3w(p.measure_f3w),
        measure_ft(p.measure_ft), measure_fw(p.measure_fw),
        measure_fl(p.measure_fl), fname_gammaw(p.fname_gammaw) {

    K = gf<imtime>{K_[0].mesh(), make_shape(n_color, n_color)};
    Kprime = gf<imtime>{K_[0].mesh(), make_shape(n_color, n_color)};
    for (auto const &t : K.mesh()) {
      for (int i = 0; i < n_color; i++) {
        auto bl_ind_i = color_to_block_and_inner_index_impl(i, gf_struct);
        for (int j = 0; j < n_color; j++) {
          auto bl_ind_j = color_to_block_and_inner_index_impl(j, gf_struct);
          int bl_ij = bl_ind_i.first * gf_struct.size() + bl_ind_j.first;
          K[t](i, j) = K_[bl_ij](t)(bl_ind_i.second, bl_ind_j.second);
          Kprime[t](i, j) = Kprime_[bl_ij](t)(bl_ind_i.second, bl_ind_j.second);
        } // j
      }   // i
    }     // t
  }
};
} // namespace triqs_ctseg
