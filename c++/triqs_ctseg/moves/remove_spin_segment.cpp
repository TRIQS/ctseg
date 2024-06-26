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

#include "remove_spin_segment.hpp"
#include "../logs.hpp"
#include <cmath>

namespace triqs_ctseg::moves {

  double remove_spin_segment::attempt() {

    LOG("\n =================== ATTEMPT REMOVE SPIN ================ \n");

    // ------------ Choose a Jperp line --------------

    auto &jl = config.Jperp_list;

    if (jl.empty()) {
      LOG("No bosonic lines!");
      return 0;
    }

    line_idx = rng(jl.size());

    // ------------- Find the hanging segment -------------

    // Find the two segments whose c is connected to the spin line
    auto it_up   = lower_bound(config.seglists[0], config.Jperp_list[line_idx].tau_Sminus);
    auto it_down = lower_bound(config.seglists[1], config.Jperp_list[line_idx].tau_Splus);

    // The hanging segment if it is in spin up
    auto spin_seg_up = segment_t{jl[line_idx].tau_Sminus, jl[line_idx].tau_Splus};

    // The hanging segment if it is in spin down
    auto spin_seg_down = segment_t{jl[line_idx].tau_Splus, jl[line_idx].tau_Sminus};

    making_full_line = *it_up == spin_seg_up and *it_down == spin_seg_down;

    if (making_full_line) {
      LOG("Making full line.");
      // Randomly choose one of the ways of making a full line
      if (rng(2) == 0) {
        spin_seg   = spin_seg_up;
        orig_color = 0;
        dest_color = 1;
        orig_it    = it_up;
        dest_it    = it_down;
        LOG("Removing segment at spin up (color 0)");
      } else {
        spin_seg   = spin_seg_down;
        orig_color = 1;
        dest_color = 0;
        orig_it    = it_down;
        dest_it    = it_up;
        LOG("Removing segment at spin down (color 1)");
      }
    } else if (*it_up == spin_seg_up) {
      spin_seg   = spin_seg_up;
      orig_color = 0;
      dest_color = 1;
      orig_it    = it_up;
      dest_it    = it_down;
      LOG("Removing segment at spin up (color 0)");
    } else if (*it_down == spin_seg_down) {
      spin_seg   = spin_seg_down;
      orig_color = 1;
      dest_color = 0;
      orig_it    = it_down;
      dest_it    = it_up;
      LOG("Removing segment at spin down (color 1)");
    } else {
      LOG("Bosonic line does not correspond to a segment.");
      return 0;
    }

    auto &dsl      = config.seglists[dest_color];
    dest_right_idx = dest_it - dsl.cbegin();
    dest_left_idx  = modulo(dest_right_idx - 1, dsl.size());

    // ------------- Check if space is free -------------

    if (dsl[dest_left_idx].tau_cdag != spin_seg.tau_c) {
      LOG("Space in opposite line is occupied.");
      return 0;
    }
    LOG("Regrouping between positions {} and {} in opposite spin.", dest_left_idx, dest_right_idx);

    // ------------  Trace ratio  -------------

    double ln_trace_ratio = (wdata.mu(dest_color) - wdata.mu(orig_color)) * spin_seg.length();
    if (wdata.has_Dt) {
      for (auto c : range(config.n_color())) {
        ln_trace_ratio -= K_overlap(config.seglists[c], spin_seg.tau_c, spin_seg.tau_cdag, wdata.K, orig_color, c);
        // "antisegment" - careful with order
        ln_trace_ratio -= K_overlap(config.seglists[c], spin_seg.tau_cdag, spin_seg.tau_c, wdata.K, dest_color, c);
      }
      // Correct for the interactions of the removed operators with themselves
      ln_trace_ratio -= real(wdata.K(double(spin_seg.length()))(orig_color, orig_color));
      ln_trace_ratio -= real(wdata.K(double(spin_seg.length()))(dest_color, dest_color));
      ln_trace_ratio += 2 * real(wdata.K(double(spin_seg.length()))(orig_color, dest_color));
    }
    double trace_ratio = std::exp(ln_trace_ratio);
    trace_ratio /= -(real(wdata.Jperp(double(spin_seg.length()))(0, 0)) / 2);

    // ------------  Det ratio  ---------------

    double det_ratio = 1.0;

    // ------------  Proposition ratio ------------

    tau_t new_seg_length = making_full_line ? tau_t::beta() : dsl[dest_left_idx].tau_c - dsl[dest_right_idx].tau_cdag;
    double future_number_seg = making_full_line ? 1 : double(dsl.size()) - 1;

    // T direct  = 1 / #Jperp * (full_line ? 1/2 : 1)  /// because if full_line, rng(2) above !
    // T inverse =  1/ (# color * #seg in future * new_seg_len^2) * (future_full_line ? 1 : 2)
    // the (splitting_full_line ...) simplify to a ratio of 2
    double prop_ratio = double(config.Jperp_list.size())
       / (double(config.n_color()) * future_number_seg * new_seg_length * new_seg_length / 2);

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double remove_spin_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    auto &sl  = config.seglists[orig_color];
    auto &dsl = config.seglists[dest_color];

    // Remove segment at origin
    sl.erase(orig_it);

    // Regroup segment at destination
    if (making_full_line) {
      dsl[0] = segment_t{tau_t::beta(), tau_t::zero()};
    } else {
      auto new_seg       = segment_t{dsl[dest_left_idx].tau_c, dsl[dest_right_idx].tau_cdag, dsl[dest_left_idx].J_c,
                               dsl[dest_right_idx].J_cdag};
      dsl[dest_left_idx] = new_seg;
      dsl.erase(dsl.begin() + dest_right_idx);
    }

    // Remove Jperp line
    auto &jl = config.Jperp_list;
    jl.erase(jl.begin() + line_idx);

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata);
    LOG("Configuration is {}", config);

    return 1.0;
  }

  //--------------------------------------------------
  void remove_spin_segment::reject() { LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n"); }

} // namespace triqs_ctseg::moves
