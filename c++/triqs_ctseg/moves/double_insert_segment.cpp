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

#include "double_insert_segment.hpp"
#include "../logs.hpp"
#include <cmath>

namespace triqs_ctseg::moves {

  double double_insert_segment::attempt() {

    LOG("\n =================== ATTEMPT DOUBLE INSERT ================ \n");

    // ------------ Choice of segment --------------
    // Select insertion colors
    int color_0 = rng(config.n_color());
    int color_1 = rng(config.n_color() - 1);
    if (color_1 >= color_0) ++color_1; // little trick to select another color
    colors.assign({color_0, color_1});
    LOG("Inserting at color ({}, {})", color_0, color_1);

    // For each color, perform a single insert
    for (auto const &[i, color] : itertools::enumerate(colors)) {
      auto &sl = config.seglists[color];

      // Select insertion window [tau_left,tau_right]
      tau_t tau_left = tau_t::beta(), tau_right = tau_t::zero();

      // If sl.empty, no segment, we have the window.
      // Otherwise, we pick up one segment at random
      if (not sl.empty()) {
        if (is_full_line(sl.back())) {
          LOG("Full line, cannot insert. (color {})", color);
          return 0;
        }
        // Randomly choose one existing segment
        long seg_idx = rng(sl.size());
        tau_left     = sl[seg_idx].tau_cdag;                     // tau_left is cdag of this segment
        tau_right    = sl[modulo(seg_idx + 1, sl.size())].tau_c; // tau_right is c of next segment, possibly cyclic
      }

      // We now have the insertion window [tau_left,tau_right]
      // Warning : it can be cyclic !
      LOG("Insertion window is tau_left = {}, tau_right = {} (color {})", tau_left, tau_right, color);
      window_length[i] = tau_left - tau_right;

      // Choose two random times in insertion window
      auto dt1 = tau_t::random(rng, window_length[i]);
      auto dt2 = tau_t::random(rng, window_length[i]);
      if (dt1 == dt2) {
        LOG("Insert_segment: generated equal times. Rejecting (color {})", color);
        return 0;
      }
      // We ensure that dt1 < dt2, unless inserting in empty line because :
      // 1- in non empty line, we insert a segment in a void so a [c cdag]
      // 2- in an empty line, we can insert [c cdag] or [cdag c]
      // BEWARE in the reverse move.
      if (dt1 > dt2 and not sl.empty()) std::swap(dt1, dt2); // if inserting into an empty line, two ways to insert

      // Here is the segment we propose to insert
      prop_seg[i] = segment_t{tau_left - dt1, tau_left - dt2};
      LOG("Inserting segment with c at {}, cdag at {} (color {})", prop_seg[i].tau_c, prop_seg[i].tau_cdag, color);
    } // color

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = 0.0;
    for (auto const &[i, color] : itertools::enumerate(colors)) {
      ln_trace_ratio += wdata.mu(color) * prop_seg[i].length(); // chemical potential
      // Overlaps
      for (auto c : range(config.n_color())) {
        if (c != color) ln_trace_ratio += -wdata.U(color, c) * overlap(config.seglists[c], prop_seg[i]);
        if (wdata.has_Dt)
          ln_trace_ratio += K_overlap(config.seglists[c], prop_seg[i].tau_c, prop_seg[i].tau_cdag, wdata.K, color, c);
      }
      if (wdata.has_Dt)
        ln_trace_ratio += -real(wdata.K(double(prop_seg[i].length()))(color, color)); // Correct double counting
    } // color

    // Counting the overlap between the inserting segments
    // Make the prop_seg[1] as a seglist to use K_overlap()
    std::vector<segment_t> seglist_temp = {prop_seg[1]};
    ln_trace_ratio += -wdata.U(color_0, color_1) * overlap(prop_seg[1], prop_seg[0]);
    if (wdata.has_Dt)
      ln_trace_ratio += K_overlap(seglist_temp, prop_seg[0].tau_c, prop_seg[0].tau_cdag, wdata.K, color_0, color_1);

    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    //  insert tau_cdag as a line (first index) and tau_c as a column (second index).
    auto &bl_0       = wdata.block_number[color_0];
    auto &bl_1       = wdata.block_number[color_1];
    auto &bl_idx_0   = wdata.index_in_block[color_0];
    auto &bl_idx_1   = wdata.index_in_block[color_1];
    is_same_block    = bl_0 == bl_1;
    double det_ratio;
    if (is_same_block) { // insert two rows and columns on the same block
      auto &D        = wdata.dets[bl_0];
      if (wdata.offdiag_Delta) {
        if (cdag_in_det(prop_seg[0].tau_cdag, D) or c_in_det(prop_seg[0].tau_c, D) or 
            cdag_in_det(prop_seg[1].tau_cdag, D) or c_in_det(prop_seg[1].tau_c, D)) {
          LOG("One of the proposed times already exists in another line of the same block. Rejecting.");
          return 0;
        }
      }
      det_ratio = D.try_insert2(det_lower_bound_x(D, prop_seg[0].tau_cdag),
                                det_lower_bound_x(D, prop_seg[1].tau_cdag),
                                det_lower_bound_y(D, prop_seg[0].tau_c),
                                det_lower_bound_y(D, prop_seg[1].tau_c),
                                {prop_seg[0].tau_cdag, bl_idx_0},
                                {prop_seg[1].tau_cdag, bl_idx_1},
                                {prop_seg[0].tau_c, bl_idx_0},
                                {prop_seg[1].tau_c, bl_idx_1});
    }
    else { // insert on different blocks
      auto &D_0 = wdata.dets[bl_0];
      auto &D_1 = wdata.dets[bl_1];
      if (wdata.offdiag_Delta) {
        if (cdag_in_det(prop_seg[0].tau_cdag, D_0) or c_in_det(prop_seg[0].tau_c, D_0) or 
            cdag_in_det(prop_seg[1].tau_cdag, D_1) or c_in_det(prop_seg[1].tau_c, D_1)) {
          LOG("One of the proposed times already exists in another line of the same block. Rejecting.");
          return 0;
        }
      }
      auto det_ratio_0 = D_0.try_insert(det_lower_bound_x(D_0, prop_seg[0].tau_cdag),
                                        det_lower_bound_y(D_0, prop_seg[0].tau_c),
                                        {prop_seg[0].tau_cdag, bl_idx_0},
                                        {prop_seg[0].tau_c, bl_idx_0});
      auto det_ratio_1 = D_1.try_insert(det_lower_bound_x(D_1, prop_seg[1].tau_cdag),
                                        det_lower_bound_y(D_1, prop_seg[1].tau_c),
                                        {prop_seg[1].tau_cdag, bl_idx_1},
                                        {prop_seg[1].tau_c, bl_idx_1});
      det_ratio        = det_ratio_0 * det_ratio_1;
    }

    // ------------  Proposition ratio ------------
    double prop_ratio = 1.0;
    for (auto const &[i, color] : itertools::enumerate(colors)) {
      auto &sl = config.seglists[color];
      double current_number_intervals = std::max(long(1), long(sl.size()));
      double future_number_segments   = sl.size() + 1;
      // T direct  = 1 / current_number_intervals * 1/window_length^2 *
      //                * (2 iif not empty as the proba to get the dt1, dt2 coupled is x 2)
      // T inverse = 1 / future_number_segments
      prop_ratio *= (current_number_intervals * window_length[i] * window_length[i] / (sl.empty() ? 1 : 2)) / future_number_segments;
    // Account for absence of time swapping when inserting into empty line.
    } // color

    LOG("trace_ratio = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;

    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double double_insert_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    double initial_sign = trace_sign(wdata);
    LOG("Initial sign is {}. Initial configuration: {}", initial_sign, config);

    // Insert the times into the det
    if (is_same_block)
      wdata.dets[wdata.block_number[colors[0]]].complete_operation();
    else
      for (auto const &color : colors)
        wdata.dets[wdata.block_number[color]].complete_operation();

    // Insert the segment in an ordered list
    for (auto const &[i, color] : itertools::enumerate(colors)) {
      auto &sl = config.seglists[color];
      sl.insert(std::upper_bound(sl.begin(), sl.end(), prop_seg[i]), prop_seg[i]);
    } // color

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata);

    double final_sign = trace_sign(wdata);
    double sign_ratio = final_sign / initial_sign;
    LOG("Final sign is {}", final_sign);

    if (sign_ratio * det_sign == -1.0) wdata.minus_sign = true;

    LOG("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void double_insert_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    if (is_same_block)
      wdata.dets[wdata.block_number[colors[0]]].reject_last_try();
    else
      for (auto const &color : colors)
        wdata.dets[wdata.block_number[color]].reject_last_try();
  }

} // namespace triqs_ctseg::moves
