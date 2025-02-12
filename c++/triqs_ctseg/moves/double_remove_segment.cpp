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

#include "double_remove_segment.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::moves {

  double double_remove_segment::attempt() {

    LOG("\n =================== ATTEMPT DOUBLE REMOVE ================ \n");

    // ------------ Choice of segment --------------
    // Select removal colors
    int color_0 = rng(config.n_color());
    int color_1 = rng(config.n_color() - 1);
    if (color_1 >= color_0) ++color_1; // little trick to select another color
    colors.assign({color_0, color_1});
    LOG("Removing at color ({}, {})", color_0, color_1);

    // For each color, perform a single remove
    for (auto const &[i, color] : itertools::enumerate(colors)) {
      auto &sl = config.seglists[color];

      // If color is empty, nothing to remove
      if (sl.empty()) {
        LOG("Empty line, cannot remove. (color {})", color);
        return 0;
      }

      // Select segment to remove
      prop_seg_idx[i] = rng(sl.size());
      prop_seg[i]     = sl[prop_seg_idx[i]];
      if (is_full_line(prop_seg[i])) {
        LOG("Cannot remove full line. (color {})", color);
        return 0;
      }
      if (prop_seg[i].J_c or prop_seg[i].J_cdag) {
        LOG("Segment has spin line attached, cannot remove. (color {})", color);
        return 0;
      }

      LOG("Removing segment at position {} : c at {}, cdag at {} (color {})", prop_seg_idx[i], prop_seg[i].tau_c, prop_seg[i].tau_cdag, color);
    } // color

    // ------------  Trace ratio  -------------
    // Same as double insert, up to the sign
    // FIXME : pull it out ?
    double ln_trace_ratio = 0.0;
    for (auto const &[i, color] : itertools::enumerate(colors)) {
      ln_trace_ratio += -wdata.mu(color) * prop_seg[i].length();
      for (auto c : range(config.n_color())) {
        if (c != color) ln_trace_ratio -= -wdata.U(color, c) * overlap(config.seglists[c], prop_seg[i]);
        if (wdata.has_Dt)
          ln_trace_ratio -= K_overlap(config.seglists[c], prop_seg[i].tau_c, prop_seg[i].tau_cdag, wdata.K, color, c);
      }
      if (wdata.has_Dt) ln_trace_ratio -= real(wdata.K(double(prop_seg[i].length()))(color, color));
    } // color

    // Counting the overlap between the removal segments
    // Make the prop_seg[1] as a seglist to use K_overlap()
    // Be careful to the sign here!
    std::vector<segment_t> seglist_temp = {prop_seg[1]};
    ln_trace_ratio += -wdata.U(color_0, color_1) * overlap(prop_seg[1], prop_seg[0]);
    if (wdata.has_Dt)
      ln_trace_ratio += K_overlap(seglist_temp, prop_seg[0].tau_c, prop_seg[0].tau_cdag, wdata.K, color_0, color_1);

    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    // same code as in double insert. In Insert, it is a true bound, does not insert at same time
    auto &bl_0       = wdata.block_number[color_0];
    auto &bl_1       = wdata.block_number[color_1];
    is_same_block = bl_0 == bl_1;
    double det_ratio;
    if (is_same_block) { // remove two rows and columns on the same block
      auto &D        = wdata.dets[bl_0];
      det_ratio = D.try_remove2(det_lower_bound_x(D, prop_seg[0].tau_cdag),
                                det_lower_bound_x(D, prop_seg[1].tau_cdag),
                                det_lower_bound_y(D, prop_seg[0].tau_c),
                                det_lower_bound_y(D, prop_seg[1].tau_c));
    }
    else { // remove on different blocks
      auto &D_0        = wdata.dets[bl_0];
      auto &D_1        = wdata.dets[bl_1];
      auto det_ratio_0 = D_0.try_remove(det_lower_bound_x(D_0, prop_seg[0].tau_cdag),
                                        det_lower_bound_y(D_0, prop_seg[0].tau_c));
      auto det_ratio_1 = D_1.try_remove(det_lower_bound_x(D_1, prop_seg[1].tau_cdag),
                                        det_lower_bound_y(D_1, prop_seg[1].tau_c));
      det_ratio        = det_ratio_0 * det_ratio_1;
    }

    // ------------  Proposition ratio ------------
    double prop_ratio = 1.0;
    for (auto const &[i, color] : itertools::enumerate(colors)) {
      auto &sl = config.seglists[color];
      double current_number_segments = sl.size();
      double future_number_intervals = std::max(1, int(sl.size()) - 1);
      // Insertion window for the reverse move double_insert_segment
      // initialise at (beta,0)
      auto tau_left = tau_t::beta(), tau_right = tau_t::zero();
      if (current_number_segments != 1) {
        // Find left, right, with cyclicity
        tau_right = sl[modulo(prop_seg_idx[i] + 1, sl.size())].tau_c;
        tau_left  = sl[modulo(prop_seg_idx[i] - 1, sl.size())].tau_cdag;
      }
      auto window_length = double(tau_left - tau_right);

      prop_ratio *= current_number_segments / 
        (future_number_intervals * window_length * window_length / (current_number_segments == 1 ? 1 : 2));
    }
    LOG("trace_ratio = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    double prod = trace_ratio * det_ratio * prop_ratio;

    return (std::isfinite(prod)) ? prod : det_sign;
  }

  //--------------------------------------------------

  double double_remove_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    double initial_sign = trace_sign(wdata);
    LOG("Initial sign is {}. Initial configuration: {}", initial_sign, config);

    // Update the dets
    if (is_same_block) 
      wdata.dets[wdata.block_number[colors[0]]].complete_operation();
    else 
      for (auto const &color : colors)
        wdata.dets[wdata.block_number[color]].complete_operation();

    // Remove the segments
    for (auto const &[i, color] : itertools::enumerate(colors)) {
      auto &sl = config.seglists[color];
      sl.erase(sl.begin() + prop_seg_idx[i]);
    }

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
  void double_remove_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    if (is_same_block)
      wdata.dets[wdata.block_number[colors[0]]].reject_last_try();
    else
      for (auto const &color : colors)
        wdata.dets[wdata.block_number[color]].reject_last_try();
  }

} // namespace triqs_ctseg::moves
