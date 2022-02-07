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
#include "segment.hpp"
#include "t_ordered_colored_c_ops.hpp"
#include "util.hpp"
#include <triqs/det_manip/det_manip.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

namespace triqs_ctseg {

/// a version which stores a new map....
struct hyb_opmap {
  typedef std::map<qmc_time_t, colored_const_iterator, std::greater<qmc_time_t>>
      ops_with_delta_lines_t;
  ops_with_delta_lines_t::const_iterator l, r;
  std::vector<ops_with_delta_lines_t> m; // one map per line/color
  gf_struct_t const &gf_struct;
  hyb_opmap(t_ordered_colored_c_ops const &ops_map,
            gf_struct_t const &gf_struct_)
      : m(ops_map.n_colors()), gf_struct(gf_struct_) {}

  void insert_in_map(colored_const_iterator const &op1,
                     colored_const_iterator const &op2) {
    int c = op1->color;
    bool ok;
    std::tie(l, ok) = m[c].insert({op1->tau, op1});
    if (!ok)
      throw insertion_error(); // CANNOT HAPPEN !
    std::tie(r, ok) = m[c].insert(std::make_pair(op2->tau, op2));
    if (!ok) {
      m[c].erase(l);
      throw insertion_error();
    }
  }

  void reject_insert_in_map(int color) {
    m[color].erase(l);
    m[color].erase(r);
  }

  void remove_from_map(colored_const_iterator const &op) {
    // std::cerr <<"---remove "<<double(op->tau)<<" from map...------" << *this
    // << std::endl;
    int c = op->color;
    // std::cerr << " col = " << c << std::endl;
    auto ll = m[c].find(op->tau);
    // std::cerr << " ll = "<< ll << std::endl;
    assert(ll != m[c].end());
    m[c].erase(ll);
  }

  void change_in_map(qmc_time_t old_tau, colored_const_iterator const &new_op) {
    std::cerr << "change in map" << *this << ", looking for " << double(old_tau)
              << std::endl;
    int c = new_op->color;
    auto ll = m[c].find(old_tau);
    assert(ll != m[c].end());
    m[c].erase(ll);

    std::cerr << "inserting..." << std::endl;
    bool ok;
    std::tie(l, ok) = m[c].insert({new_op->tau, new_op});
    if (!ok)
      throw insertion_error(); // CANNOT HAPPEN !
  }
  /// return the index of the operator, counting only the operators of the same
  /// kind and same block_index before [creation|annihilation]
  size_t find_index(colored_const_iterator const &it) const {
    int block_index =
        color_to_block_and_inner_index_impl(it->color, gf_struct).first;
    size_t i = 0;
    // count #ops of same kind and same block *before* it
    for (int c = 0; c < it->color; c++)
      if (color_to_block_and_inner_index_impl(c, gf_struct).first ==
          block_index)
        i += seg_number(c);
    for (auto pos = begin(m[it->color]); pos->first != it->tau; ++pos)
      if (pos->second->dagger == it->dagger)
        ++i;
    return i;
  }

  /// number of hybridization lines /segments on line 'color'
  int seg_number(int color) const { return m[color].size() / 2; }

  /// given a color and a int, returns the corresponding *hybridized* segment
  segment find(int color, const size_t op_index) const {
    auto l = begin(m[color]);
    for (int i = 0; i < op_index; ++i) {
      ++l;
    }
    if (l == end(m[color]))
      TRIQS_RUNTIME_ERROR << "segment of index " << op_index << " on line "
                          << color << " not found";
    auto r = l;
    ++r;
    if (r == end(m[color]))
      r = begin(m[color]); // a cyclic right
    return {l->second, r->second};
  }

  int antiseg_number(int color) const {
    auto l = begin(m[color]);
    int n = 0;
    for (int i = 0; i < seg_number(color); ++i) {
      if (l->second->dagger)
        n++;
      ++l;
      ++l;
    }
    return n;
  }

  friend std::ostream &operator<<(std::ostream &out, hyb_opmap const &h) {
    int i = 0;
    auto const &m = h.m;
    for (auto it = m.begin(); it != m.end(); it++, i++) {
      out << "[flav #" << i << ":";
      for (auto it2 = it->begin(); it2 != it->end(); it2++)
        out << double(it2->first) << " - ";
      out << "] ";
    }
    return out;
  }
  friend void h5_write(h5::group fg, std::string subgroup_name,
                       hyb_opmap const &h) {
    h5::group gr = fg.create_group(subgroup_name);
    auto const &m = h.m;
    int tot_number = 0;
    for (auto it = m.begin(); it != m.end(); it++)
      for (auto it2 = it->begin(); it2 != it->end(); it2++)
        tot_number++;
    if (tot_number > 0) {
      array<double, 2> a(tot_number, 2);
      int k = 0, i = 0;
      for (auto it = m.begin(); it != m.end(); it++, i++) {
        for (auto it2 = it->begin(); it2 != it->end(); it2++, k++) {
          a(k, 0) = double(it2->first); // time
          a(k, 1) = i;                  // color
        }
      }
      h5_write(gr, "hybtime_table", a);
    }
  }
};

} // namespace triqs_ctseg
