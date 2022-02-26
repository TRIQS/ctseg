#include "work_data.hpp"
#include <triqs/operators/util/extractors.hpp>

work_data_t::work_data_t(param_t const &p, inputs_t const &inputs) {

  beta = p.beta;

  // ......... U .................
  n_color = 0;
  // std::vector<std::string> block_names;
  for (auto const &[bl_name, bl_size] : p.gf_struct) {
    // block_names.push_back(bl_name);
    n_color += bl_size;
  }

  // color dependent chemical potential
  nda::vector<double> mu(n_color);
  mu = 0;

  if (p.hartree_shift.size() > 0) {
    ALWAYS_EXPECTS((p.hartree_shift.size() == n_color), "Hartree shift size is not {}", n_color);
    mu += p.hartree_shift;
  }

  // ......... U .................
  // Extract the U from the operator
  auto U_full           = triqs::operators::utils::dict_to_matrix(triqs::operators::utils::extract_U_dict2(p.h_int), gf_struct);
  nda::matrix<double> U = real(U_full);

  // report
  if (c.rank() == 0) {
    spdlog::info("mu = {}\n U = {}", mu, U);
    spdlog::info("dynamical_U = {}\n jperp_interactions = {}\n ", has_Dt, has_jperp);
  }

  // Check
  ALWAYS_EXPECTS((jperp_interactions and deltaw.size() != 2), "Error : jperp interactions is true and we have {} blocks instead of 2", deltaw.size());

  // .........  Do we have dynamical interactions ?
  has_Dt    = !is_zero(d0t);
  has_jperp = !is_zero(jperpt);

  // FIXME : K, Kprime .. + adjust.

  // ........   dets .....................

  delta = map([](gf_const_view<imtime> d) { return real(d); }, inputs.delta);

  for (auto const &bl : range(delta.size())) {
    // FIXME : spdlog
    if (!is_gf_real(delta[bl], 1e-10)) {
      std::cerr << "WARNING: The Delta(tau) block number " << bl << " is not real in tau space\n";
      std::cerr << "WARNING: max(Im[Delta(tau)]) = " << max_element(abs(imag(delta[bl].data()))) << "\n";
      std::cerr << "WARNING: Dissregarding the imaginary component in the calculation.\n";
    }
    dets.emplace_back(delta_block_adaptor(real(delta[bl])), p.det_init_size);
    dets.back().set_singular_threshold(p.det_singular_threshold);
    dets.back().set_n_operations_before_check(p.det_n_operations_before_check);
    dets.back().set_precision_warning(p.det_precision_warning);
    dets.back().set_precision_error(p.det_precision_error);
  }
}


#if 0 
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
	}  j
      }    i
    }      t */

  if (dynamical_U) {
    // extract kt and kprimet from d0w
    for (int bl = 0; bl < kt.size(); bl++) {
      for (auto const &t : kt[bl].mesh()) {
        for (auto i = 0; i < kt[bl].target_shape()[0]; i++) {
          for (auto j = 0; j < kt[bl].target_shape()[1]; j++) {
            kt[bl][t](i, j) = 0.0;
            kprimet[bl][t](i, j) = 0.0;

            for (auto const &w : d0w[bl].mesh()) {
              if (w.n > 0) {
                double mats = dcomplex(w).imag();
                kt[bl][t](i, j) += real(d0w[bl](w)(i, j)) / (mats * mats) *
                                   (1 - cos(mats * t));

                kprimet[bl][t](i, j) +=
                    real(d0w[bl](w)(i, j)) / mats * sin(mats * t);
              }
            }
            kt[bl][t](i, j) *= 2. / beta;
            kprimet[bl][t](i, j) *= 2. / beta;
            kt[bl][t](i, j) +=
                0.5 * real(d0w[bl](0)(i, j)) * t * (t / beta - 1.);
            kprimet[bl][t](i, j) += real(d0w[bl](0)(i, j)) * (t / beta - .5);
          }
        }
      }
    }
    // extract kt and kprimet from d0w
    for (auto &t : kperpprimet.mesh()) {
      kperpprimet[t] = 0.0;
      for (auto &w : jperpw.mesh()) {
        double mats = dcomplex(w).imag();
        if (w.n > 0)
          kperpprimet[t] += real(jperpw(w)) / mats * sin(mats * t);
      }
      kperpprimet[t] *= 2. / beta;
      kperpprimet[t] += real(jperpw(0)) * (t / beta - .5);
      kperpprimet[t] *= 0.25; // D^z = J_z/4.; K pertains to D, not J
    }

    // FIXME CHECK 
    gf<imfreq> D_sp_minus_one_fourth_Jperp(jperpw.mesh(), {1, 1});
    auto ind_UV = [&](std::string const &s) {
      return get_index(d0w.block_names(), s);
    };
    for (auto const &iom : D_sp_minus_one_fourth_Jperp.mesh())
      D_sp_minus_one_fourth_Jperp[iom] =
          0.5 *
              (d0w[ind_UV(block_names[0] + "|" + block_names[0])](iom)(0, 0) -
               d0w[ind_UV(block_names[0] + "|" + block_names[1])](iom)(0, 0)) -
          0.25 * jperpw(iom)(0, 0);
    full_spin_rot_inv =
        is_zero(D_sp_minus_one_fourth_Jperp) and jperp_interactions;
  } // end if dynamical_U

  // KEEP 
  // Get the new values for mu and Umatrix from Kprime
  if (dynamical_U) {
    for (size_t i = 0; i < U.shape()[0]; ++i) {
      auto bl_ind_i = color_to_block_and_inner_index_impl(i, gf_struct);
      auto bl_ii = bl_ind_i.first * gf_struct.size() +
                   bl_ind_i.first; // double index (e.g. up|up)
      mu(i) +=
          kprimet[bl_ii][closest_mesh_pt(0.0)](bl_ind_i.second, bl_ind_i.second)
              .real();
      for (size_t j = 0; j < U.shape()[1]; ++j) {
        auto bl_ind_j = color_to_block_and_inner_index_impl(j, gf_struct);
        auto bl_ij = bl_ind_i.first * gf_struct.size() + bl_ind_j.first;
        if (i != j)
          U(i, j) -= 2. * kprimet[bl_ij][closest_mesh_pt(0.0)](bl_ind_i.second,
                                                               bl_ind_j.second)
                              .real();
      }
    }
  }
#endif
