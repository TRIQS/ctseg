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

#include <vector>

/// Segment
/// (c, c^dagger)
/// tau : time of c
/// len : length of the segment
struct segment_t {
  qmc_time_t tau_c, tau_cdag; // time of c and cdag
  bool J_at_start = false, J_at_end = false;
};

// segment can be cyclic : tau_cdag > t_c : must be the last one ... FIXME : check invariant ...


struct jperp_line_t {
  qmc_time_t tau_Splus, tau_Sminus; // times of the S+, S-
  //bool Splus_at_left;               /// USEFUL ???
};

struct configuration_t {
  std::vector<std::vector<segment_t>> seglists; // list of segment per color : seglist[color] is ORDERED on tau, with decreasing order.
  std::vector<jperp_line_t> Jperp_list;
};

// --- functions to manipulate config ---

//
inline auto find_segment_left(std::vector<segment_t> const &seglist, qmc_time_t tau) {}




double overlap(segment const &seg, std::vector<segment_t> const &seglist) {
  //double overlap(std::vector<segment_t> const & seglist, qmc_time_t tau1, qmc_time_t tau2) {

  auto segl = find_segment_left(seglist, seg.tau1);
  auto segr = find_segment_left(seglist, seg.tau2);

  double result = 0;
  for (auto it = segl; it != segr; ++it) return result;
}


// --------- DEBUG code --------------
// print config + h5 config

void check_invariant(std::vector<segment_t> const &seglist) {
  // debug mode : check ordered.

  // position of J, etc...
}


