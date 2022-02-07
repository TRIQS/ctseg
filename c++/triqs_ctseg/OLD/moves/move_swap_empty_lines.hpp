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
#include "../segment.hpp"

namespace triqs_ctseg {
/// Swap empty lines
/**
 *Move which tries swapping an empty line with a full line if no operator on
 *line.
 */
class move_swap_empty_lines {
  qmc_parameters *params;
  configuration *config;
  triqs::mc_tools::random_generator &RND;
  int c1, c2;
  bool tried_swap;

public:
  move_swap_empty_lines(qmc_parameters *params_, configuration *config_,
                        triqs::mc_tools::random_generator &RND_)
      : params(params_), config(config_), RND(RND_){};

  double attempt() {
#ifdef CTSEG_DEBUG
    if (config->print_condition())
      std::cout << "\n =================== ATTEMPT SWAP EMPTY LINES "
                   "================ \n"
                << std::endl;
#endif
    tried_swap = false;
    // choose two lines: one original color, one final color
    c1 = RND(params->n_color);
    c2 = RND(params->n_color - 1);

    // reject in most cases
    if (c2 == c1)
      return 0.0;
    if ((config->ops_map().seg_number(c1) > 0) ||
        (config->ops_map().seg_number(c2) > 0))
      return 0.0;
    if (config->trace.full_lines(c1) == config->trace.full_lines(c2))
      return 0.0; // the lines are both full or both empty

    double ln_trace_ratio = config->trace.swap_full_empty(c1, c2);
    tried_swap = true;
    return std::exp(ln_trace_ratio); // det_ratio and prop_ratio are 1
  }

  //-----------------------------------

  double accept() {
    config->id++;
#ifdef CTSEG_DEBUG
    if (config->print_condition())
      std::cout << "- - - - - ====> ACCEPT - - - - - - - - - - -" << std::endl;
    std::cerr << "config: " << *config << std::endl;
#endif
    return 1.0;
  }

  //-----------------------------------

  void reject() {
    config->id++;
    if (tried_swap)
      config->trace.swap_full_empty(c1, c2);
#ifdef CTSEG_DEBUG
    if (config->print_condition())
      std::cout << "- - - - - ====> REJECT - - - - - - - - - - -" << std::endl;
    std::cerr << "config: " << *config << std::endl;
#endif
  }
};

} // namespace triqs_ctseg
