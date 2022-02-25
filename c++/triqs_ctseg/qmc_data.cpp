#include "qmc_data.hpp"

 // Construction
qmc_data::qmc_data(params_t const &p, 
    block_gf_const_view<imtime> delta delta(map([](gf_const_view<imtime> d) { return real(d); }, delta // FIXME!!!
    gf<imtime> const & d0t, gf<imtime> const & jperpt_)
    {

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
    

/*     K      = gf<imtime>{K_[0].mesh(), make_shape(n_color, n_color)};
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
    }     // t */
  }