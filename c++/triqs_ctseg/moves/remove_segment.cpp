// Copyright (c) 2022-2024 Simons Foundation
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Nikita Kavokine, Hao Lu, Olivier Parcollet, Nils Wentzell

#include "remove_segment.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::moves {

  double remove_segment::attempt() {

    LOG("\n =================== ATTEMPT REMOVE ================ \n");

    // ------------ Choice of segment --------------
    // Select removal color
    color    = rng(config.n_color());
    auto &sl = config.seglists[color];
    LOG("Removing at color {}", color);

    // If color is empty, nothing to remove
    if (sl.empty()) {
      LOG("remove_segment: reject : color is empty.");
      return 0;
    }

    // Select segment to remove
    prop_seg_idx = rng(sl.size());
    prop_seg     = sl[prop_seg_idx];
    if (is_full_line(prop_seg)) {
      LOG("Cannot remove full line.");
      return 0;
    }
    if (prop_seg.J_c or prop_seg.J_cdag) {
      LOG("Segment has spin line attached, cannot remove.");
      return 0;
    }

    LOG("Removing segment at position {} : c at {}, cdag at {}", prop_seg_idx, prop_seg.tau_c, prop_seg.tau_cdag);

    // ------------  Trace ratio  -------------
    // Same as insert, up to the sign
    // FIXME : pull it out ?
    double ln_trace_ratio = -wdata.mu(color) * prop_seg.length();
    for (auto c : range(config.n_color())) {
      if (c != color) { ln_trace_ratio -= -wdata.U(color, c) * overlap(config.seglists[c], prop_seg); }
      if (wdata.has_Dt)
        ln_trace_ratio -= K_overlap(config.seglists[c], prop_seg.tau_c, prop_seg.tau_cdag, wdata.K, color, c);
    }
    if (wdata.has_Dt) ln_trace_ratio -= real(wdata.K(double(prop_seg.length()))(color, color));

    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    // same code as in insert. In Insert, it is a true bound, does not insert at same time
    auto bl        = wdata.block_number[color];
    auto &D        = wdata.dets[bl];
    auto det_ratio = D.try_remove(det_lower_bound_x(D, prop_seg.tau_cdag), //
                                  det_lower_bound_y(D, prop_seg.tau_c));

    // ------------  Proposition ratio ------------

    double current_number_segments = sl.size();
    double future_number_intervals = std::max(1, int(sl.size()) - 1);
    // Insertion window for the reverse move insert_segment
    // initialise at (beta,0)
    auto tau_left = tau_t::beta(), tau_right = tau_t::zero();
    if (current_number_segments != 1) {
      // Find left, right, with cyclicity
      tau_right = sl[modulo(prop_seg_idx + 1, sl.size())].tau_c;
      tau_left  = sl[modulo(prop_seg_idx - 1, sl.size())].tau_cdag;
    }
    auto window_length = double(tau_left - tau_right);

    double prop_ratio = current_number_segments
       / (future_number_intervals * window_length * window_length / (current_number_segments == 1 ? 1 : 2));
    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    double prod = trace_ratio * det_ratio * prop_ratio;

    return (std::isfinite(prod)) ? prod : det_sign;
  }

  //--------------------------------------------------

  double remove_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    double initial_sign = trace_sign(wdata);
    LOG("Initial sign is {}. Initial configuration: {}", initial_sign, config);

    // Update the dets
    wdata.dets[wdata.block_number[color]].complete_operation();

    auto &sl = config.seglists[color];
    // Remove the segment
    sl.erase(sl.begin() + prop_seg_idx);

    double final_sign = trace_sign(wdata);
    double sign_ratio = initial_sign / final_sign;
    LOG("Final sign is {}", final_sign);

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata);

    if (sign_ratio * det_sign == -1.0) wdata.minus_sign = true;

    LOG("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void remove_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[wdata.block_number[color]].reject_last_try();
  }

} // namespace triqs_ctseg::moves
