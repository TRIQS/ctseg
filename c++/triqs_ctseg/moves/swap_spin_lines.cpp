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

#include "swap_spin_lines.hpp"
#include "../logs.hpp"
#include <cmath>

namespace triqs_ctseg::moves {

  double swap_spin_lines::attempt() {

    LOG("\n =================== ATTEMPT SWAP SPIN LINES ================ \n");

    // ------------ Choose two spin lines --------------

    auto &jl = config.Jperp_list;

    if (jl.empty() or jl.size() == 1) {
      LOG("Nothing to swap!");
      return 0;
    }

    first_line_idx  = rng(jl.size());
    second_line_idx = rng(jl.size() - 1);
    if (second_line_idx >= first_line_idx) ++second_line_idx; // trick to select second_line_idx ! = first_line_idx

    auto &l1 = jl[first_line_idx];
    auto &l2 = jl[second_line_idx];

    // ------------  Trace ratio  -------------

    double J_current =                                                 //
       real(wdata.Jperp(double(l1.tau_Sminus - l1.tau_Splus))(0, 0)) * //
       real(wdata.Jperp(double(l2.tau_Sminus - l2.tau_Splus))(0, 0));

    double J_future =                                                  //
       real(wdata.Jperp(double(l1.tau_Sminus - l2.tau_Splus))(0, 0)) * //
       real(wdata.Jperp(double(l2.tau_Sminus - l1.tau_Splus))(0, 0));

    double trace_ratio = J_future / J_current;

    // ------------  Det ratio  ---------------

    double det_ratio = 1.0;

    // ------------  Proposition ratio ------------

    double prop_ratio = 1.0;

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    return (std::isfinite(prod) ? prod : 1);
  }

  // --------------------------------------------

  double swap_spin_lines::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    auto &jl = config.Jperp_list;
    auto &l1 = jl[first_line_idx];
    auto &l2 = jl[second_line_idx];

    // Swap lines
    std::swap(l1.tau_Splus, l2.tau_Splus);

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata);
    LOG("Configuration is {}", config);

    return 1.0;
  }

  //--------------------------------------------------
  void swap_spin_lines::reject() { LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n"); }

} // namespace triqs_ctseg::moves
