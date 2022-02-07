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
#include "./t_ordered_colored_c_ops.hpp"
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

namespace triqs_ctseg {
using triqs::gfs::closest_mesh_pt;
using colored_const_iterator = t_ordered_colored_c_ops::colored_const_iterator;

// returns true if tau_1 and tau_2 are on the same side with respect to tau_R
inline bool on_same_side(qmc_time_t const &tau_1, qmc_time_t const &tau_2,
                         qmc_time_t const &tau_R) {
  return (tau_1 > tau_R) == (tau_2 > tau_R);
}

// description of a segment : 2 times and a color ...
// this info is sufficient to insert it
struct segment_desc {
  qmc_time_t tl, tr;
  int color;
};

// a segment as inserted in the colored list :
// two iterators l, r to left|right operator of the segment
// l and r MUST be neighbours.
struct segment {
  colored_const_iterator l, r;
  // returns the length of the segment. takes care of beta-periodifity. <0 for
  // antiseg
  double length() const {
    double le = double(l->tau - r->tau);
    return (l->dagger ? -le : le);
  }
  colored_const_iterator a() const { return l->dagger ? r : l; } // annihilation
  colored_const_iterator c() const { return l->dagger ? l : r; } // creation
};

} // namespace triqs_ctseg
