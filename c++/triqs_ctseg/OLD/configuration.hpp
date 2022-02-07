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

#include "./bosonic_lines.hpp"
#include "./hybridization_dets.hpp"
#include "./trace_c_ops.hpp"
#include "./util.hpp"

namespace triqs_ctseg {

/// A MC configuration
struct configuration {
  int id;
  trace_c_ops trace;           // handles the trace (and the list of operators)
  hybridization_dets hyb_dets; // handles the hybridization lines
  bosonic_lines boson_lines;   // handles the bosonic lines
  configuration(const qmc_parameters *params, gf_struct_t const &gf_struct__,
                block_gf<imtime> const &hyb_raw, gf<imtime> const &Jperp,
                bool keep_Jperp_negative)
      : id(0), trace(params, &id),
        hyb_dets(trace.ops_map(), gf_struct__, hyb_raw),
        boson_lines(Jperp, keep_Jperp_negative) {
    trace.hyb_dets = &hyb_dets;
  };

  /// the time ordered list of colored c, cdag operators
  t_ordered_colored_c_ops const &ops_map() const { return trace.ops_map(); }

  gf_struct_t const &gf_struct() const { return hyb_dets.gf_struct(); }
  int n_blocks() const { return gf_struct().size(); }
  int n_colors() const { return ops_map().n_colors(); }
  bool print_condition() const { return true; }
  auto det(int k) const DECL_AND_RETURN(hyb_dets.det(k));

  std::pair<int, int>
  color_to_block_and_inner_index(int const color_number) const {
    return color_to_block_and_inner_index_impl(color_number, gf_struct());
  }

  int block_and_inner_index_to_color(int block, int inner_index) const {
    return block_and_inner_index_to_color_impl(block, inner_index, gf_struct());
  };

  int edge_state(int color) const { // returns the state (0 or 1) at time tau=0
    return trace.ops_map().op_number(color) == 0
               ? trace.full_lines(color)
               : trace.ops_map()
                     .begin(color)
                     ->dagger; // implicit cast bool to int
  }

  friend std::ostream &operator<<(std::ostream &out, configuration const &c) {
    mpi::communicator w;
    out << "[" << w.rank() << "] #" << c.id << ":" << std::endl;
    out << c.trace << std::endl;
    out << "HYB LINES: " << c.hyb_dets << std::endl;
    out << "BOSON LINES: " << c.boson_lines << std::endl;
    return out;
  }

  friend void h5_write(h5::group fg, std::string subgroup_name,
                       configuration const &c) {
    h5::group gr = fg.create_group(subgroup_name);

    h5_write(gr, "trace", c.trace);
    h5_write(gr, "hyb_dets", c.hyb_dets);
    h5_write(gr, "boson_lines", c.boson_lines);
    h5_write(gr, "n_lines", c.n_colors());
  }

  void print_to_h5() {
    std::string filename = "configs.h5";
    h5::file hfile(filename.c_str(), 'a');
    h5_write(hfile, "c_" + str(this->id), *this);
    hfile.close();
  }
};
} // namespace triqs_ctseg
