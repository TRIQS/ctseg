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
#include "./segment.hpp"
#include "./t_ordered_colored_c_ops.hpp"
#include <triqs/det_manip/det_manip.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

namespace triqs_ctseg {

// a simple version, which does NOT store a new map for
struct hyb_opmap_no_storage {
  t_ordered_colored_c_ops &m;

  hyb_opmap_no_storage(t_ordered_colored_c_ops &ops_map) : m(ops_map) {}

  void insert_in_map(colored_const_iterator const &op1,
                     colored_const_iterator const &op2) {}
  void reject_insert_in_map(int color) {}
  void remove_from_map(colored_const_iterator const &op) {}

  /// return the index of the operator on a given line, counting only the
  /// operators of the same kind [creation|annihilation]
  size_t find_index(colored_const_iterator const &it) const {
    size_t i = 0;
    for (auto pos = m.begin(it->color); pos->tau != it->tau; ++pos) {
      if (pos->dagger == it->dagger)
        ++i;
    }
    return i;
  }

  /// number of hybridization lines /segments on line 'color'
  int seg_number(int color) const { return m.seg_number(color); }

  /// given a color and a int, returns the corresponding *hybridized* segment
  segment find(int color, const size_t op_index) const {
    auto l = m.begin(color);
    for (int i = 0; i < op_index; ++i) {
      ++l;
    }
    if (l == m.end(color))
      TRIQS_RUNTIME_ERROR << "segment of index " << op_index << " on line "
                          << color << " not found";
    auto r = m.cyclic_right(l);
    return {l, r};
  }

  friend std::ostream &operator<<(std::ostream &out,
                                  hyb_opmap_no_storage const &h) {
    int i = 0;
    auto const &m = h.m;
    for (auto it = m.begin(); it != m.end(); it++, i++) {
      out << "[flav #" << i << ":";
      out << "] ";
    }
    return out;
  }
};
} // namespace triqs_ctseg
