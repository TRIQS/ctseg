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
#include "../qmc_parameters.hpp"

namespace triqs_ctseg {
/// Insert a segment
/**
 *Insert a segment $c(\tau)c^\dagger(\tau')$
 */
class move_insert_segment {
  const qmc_parameters *params;
  configuration *config;
  triqs::mc_tools::random_generator &RND;
  int color;

public:
  move_insert_segment(const qmc_parameters *params_, configuration *config_,
                      triqs::mc_tools::random_generator &RND_)
      : params(params_), config(config_), RND(RND_){};

  double attempt() {
#ifdef CTSEG_DEBUG
    if (config->print_condition())
      std::cout << "\n =================== ATTEMPT INSERT ================ \n"
                << std::endl;
#endif
    // --------  Selection of the times and color for the segment insertion
    // ------------

    // pick up color and first time
    color = RND(params->n_color);
    auto tau_1 =
        params->tau_seg.get_random_pt(RND); // choose tau between 0 and beta
    qmc_time_t lmax;
    bool need_swap;
    qmc_time_t tau_2;

    // pick up second time:
    if (config->ops_map().op_number(color) == 0) { // no operators
      tau_2 = params->tau_seg.get_random_pt(RND);
      need_swap = !(config->trace.full_lines(
          color)); // by convention, if empty/full line: tau_1 is creation,
                   // tau_2 is annihilation
    } else {
      auto R = config->ops_map().right_neighbor(
          tau_1, color); // point at right of tau_1
      assert(config->ops_map().cyclic_left(R)->tau != R->tau);
      lmax = (config->ops_map().cyclic_left(R)->tau - R->tau);
      tau_2 = params->tau_seg.get_random_pt(RND, lmax) + R->tau;
      need_swap = on_same_side(tau_1, tau_2, R->tau) ? !(tau_1 > tau_2)
                                                     : (tau_1 > tau_2);
    }
    if (need_swap)
      std::swap(tau_1, tau_2);
    auto sd = segment_desc{tau_1, tau_2, color};

    // ----    compute ratio of trace and det  ------------
#ifdef CTSEG_DEBUG
    if (config->print_condition()) {
      std::cerr << "tau_2 = " << tau_2 << std::endl;
    }
#endif
    segment seg;
    double ln_trace_ratio;
    try {
      std::tie(ln_trace_ratio, seg) = config->trace.try_insert_segment(sd);
    } catch (insertion_error const &e) {
      std::cerr << "insert error -- move insert" << std::endl;
      return 0;
    }

    double trace_ratio = std::exp(ln_trace_ratio);
    double det_ratio = config->hyb_dets.try_add(seg.l, seg.r);
    int n_ops_after = 2 * config->hyb_dets.seg_number(
                              color); // try_insert updates hyb op map size
    // ----    compute proposition ratio ------------
    double prop_ratio =
        config->ops_map().seg_number(color) == 1
            ? params->beta * params->beta / 2.0
            : params->beta * lmax /
                  (2 * n_ops_after); // compute proposal probability ratio

#ifdef CTSEG_DEBUG
    if (config->print_condition()) {
      if (seg.l->dagger) {
        std::cerr << "Trying to insert segment tau_a = " << seg.r->tau
                  << " tau_c = " << seg.l->tau << " on line " << color
                  << std::endl;
      } else {
        std::cerr << "Trying to insert segment tau_a = " << seg.l->tau
                  << " tau_c = " << seg.r->tau << " on line " << color
                  << std::endl;
      }
      std::cerr << "***RATIOS: TR= " << trace_ratio << "|DET= " << det_ratio
                << "|PROP = " << prop_ratio << std::endl;
    }
#endif

    double s = (det_ratio > 0.0 ? 1.0 : -1.0);
    double prod = trace_ratio * det_ratio * prop_ratio;
    return (std::isfinite(prod) ? prod : s);
  }

  //--------------------------------------------------
  double accept() {
#ifdef CTSEG_DEBUG
    if (config->print_condition())
      std::cout << "- - - - - ====> ACCEPT - - - - - - - - - - -" << std::endl;
#endif
    config->id++;
    config->hyb_dets.complete_add();
    double sign_ratio = config->trace.complete_insert_segment();
#ifdef CTSEG_DEBUG
    if (config->print_condition()) {
      // config->print_to_h5();
      std::cerr << "config: " << *config << std::endl;
    }
    config->trace.check_overlap_matrix_from_scratch();
#endif
    return sign_ratio;
  }

  //--------------------------------------------------
  void reject() {
#ifdef CTSEG_DEBUG
    if (config->print_condition())
      std::cout << "- - - - - ====> REJECT - - - - - - - - - - -" << std::endl;
#endif
    config->id++;
    config->hyb_dets.reject_add();
    config->trace.reject_insert_segment();
#ifdef CTSEG_DEBUG
    if (config->print_condition()) {
      // config->print_to_h5();
      std::cerr << "config: " << *config << std::endl;
    }
    config->trace.check_overlap_matrix_from_scratch();
#endif
  }
};
} // namespace triqs_ctseg
