/*******************************************************************************
 *
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
 *
 ******************************************************************************/
#pragma once
#include "./hybridization_dets.hpp"
#include "./qmc_parameters.hpp"
#include "./segment.hpp"
#include "./t_ordered_colored_c_ops.hpp"
#include <nda/blas/dot.hpp>

namespace triqs_ctseg {
using triqs::gfs::closest_mesh_pt;

/*
 * This object handles the trace.
 * It contains the list of colored fundamental operators, with a const accessor
 * to it (ops_map). It provides several methods to modify the trace :
 *  - insert_segment : to insert a segment or an antisegment (i.e. a c^+ c or a
 * c c^+).
 *  - remove_segment : to remove them
 *  - swap_full_empty : to swap two lines without operators...
 *
 *  Each operation (insert/remove) can be decomposed in a try/complete/reject,
 *  or called directly (doing a try-complete at once).
 *  The trace object has a state (neutral, tried_insert, tried_remove) to check
 * it is used properly.
 *
 *  The operations returns the log of the trace ratio between the configuration
 * before and after the operation together with the ratio of signs of the trace
 * and a handle on the inserted/removed segment to make inverse operation easy.
 *
 *  When lines have no operators, the bool full_lines(color) return whether this
 * is a full or an empty line. This has no meaning when there is any operator
 * with the given color. The state of empty lines is changed by swap_full_empty.
 *
 *  The object also computes the overlap matrix between segment (with a const
 * accessor overlap_matrix()).
 */
class trace_c_ops {

  const qmc_parameters *params;
  int *id;
  t_ordered_colored_c_ops
      ops; // the time ordered list of colored c, cdag operators
  std::vector<bool> _full_lines; // full_lines(color) = true if line 'color' is
                                 // full    // was private before...
  matrix<double>
      _overlap_matrix;   // matrix of the overlaps;diagonal terms:lengths
  std::vector<int> sign; // sign of each line

  // working data kept between try, insert
  vector<double> overlaps; // overlap of one segment with all other colors
  int color;
  segment seg;
  bool is_segment;
  double length;
  segment_desc s_desc;

  enum class status_t {
    neutral,
    tried_insert,
    tried_remove
  }; // state of the object
  status_t status;

public:
  trace_c_ops(const qmc_parameters *params_, int *id);

  const hybridization_dets *hyb_dets;
  double const &overlap_matrix(int i, int j) const {
    return _overlap_matrix(i, j);
  }
  bool full_lines(int color) const { return _full_lines[color]; }
  t_ordered_colored_c_ops const &ops_map() const { return ops; }

  //***************   Segment insertion ******************************

  /*
   * Insert a segment
   * s : description of the segment to insert
   * Returns (ln_trace_ratio, sign_ratio, segment)
   */
  std::tuple<double, double, segment> insert_segment(segment_desc const &s) {
    auto res = try_insert_segment(s);
    auto sign_ratio = complete_insert_segment();
    return std::make_tuple(res.first, sign_ratio, res.second);
  }

  /*
   * Insert a segment [try part]
   * s : description of the segment to insert
   * Returns (ln_trace_ratio, segment)
   * NB : This DOES insert the operator in the colored operator list
   */
  std::pair<double, segment> try_insert_segment(segment_desc const &s);

  /*
   * Insert a segment [complete part]
   * Returns sign_ratio
   */
  double complete_insert_segment();

  /*
   * Cancel the current insertion of a segment.
   * NB : can be also called from neutral state, i.e. when no try has been
   * called before.
   */
  void reject_insert_segment() {
    if (status == status_t::neutral)
      return;
    assert(status == status_t::tried_insert);
    rm_segment_and_update(); // erase operators from maps
    status = status_t::neutral;
  }

  //***************  Segment removal ******************************

  /*
   * Remove a segment s
   * Returns (ln_trace_ratio, sign_ratio, segment_description)
   * NB : segment_description to be able to reinsert the segment ...
   */
  std::tuple<double, double, segment_desc> remove_segment(segment const &s) {
    auto ratio = try_remove_segment(s);
    auto r = complete_remove_segment();
    return std::make_tuple(ratio, r.first, r.second);
  }

  /*
   * Remove a segment s [try part]
   * Returns ln_trace_ratio
   * NB : does NOT remove the operator from the list ....
   */
  double try_remove_segment(segment const &s);

  /*
   * Remove a segment s [complete part]
   * Returns (sign_ratio, segment_desc)
   */
  std::pair<double, segment_desc> complete_remove_segment();

  /*
   * Cancel the segment removal
   * NB : can be also called from neutral state, i.e. when no try has been
   * called before.
   */
  void reject_remove_segment() {
    if (status == status_t::neutral)
      return;
    assert(status == status_t::tried_remove);
    // nothing to do
    status = status_t::neutral;
  }

  //**************Operator shift *******************

  /// shift operator c_old to tau_new
  /**
   *
   * @return pair{trace ratio,sign_ratio, new operator}
   */
  // std::tuple<double, double, colored_const_iterator>
  // shift_operator(colored_const_iterator const & c_old, qmc_time_t tau_new);
  std::pair<double, colored_const_iterator>
  try_shift_operator(colored_const_iterator const &c_old, qmc_time_t tau_new);

  double complete_shift_operator();
  /// shift operator c_old to tau_new
  /**
   * first removes seg (c_old, R) and then inserts seg (c_new, R)
   * @return pair{trace ratio, sign ratio, new operator}
   * @note does not work because removing nonadjacent ops forbidden
   */
  /*
  std::tuple<double, double, colored_const_iterator>
  shift_operator(colored_const_iterator const & c_old, qmc_time_t tau_new){ auto
  R_old = ops.cyclic_right(c_old); int c = c_old->color;

   double ln_trace_ratio_rm, sign_ratio_rm;
   segment_desc old_seg_desc;
   std::tie(ln_trace_ratio_rm, sign_ratio_rm, old_seg_desc) =
  remove_segment(segment{c_old, R_old});

   auto new_seg_desc = segment_desc{tau_new, old_seg_desc.tr, c};

   double ln_trace_ratio_ins, sign_ratio_ins;segment new_seg;
   std::tie(ln_trace_ratio_ins, sign_ratio_ins, new_seg) =
  insert_segment(new_seg_desc);

   return std::make_tuple(ln_trace_ratio_rm + ln_trace_ratio_ins,
  sign_ratio_rm*sign_ratio_ins,  new_seg.l);
  }*/

  //*********************************************

  /*
   * Change the state of line without operators  full <-> empty for color c
   * Returns : ln_trace_ratio
   * Precondition : there must not be any operator of color c, or it is a
   * logical error in the code.
   */
  double swap_full_empty(int c1, int c2);

  void rm_segment_and_update();

  // sign of a line: flip sign if odd number of segments and starts with dagger
  //@note supposes hyb_dets are correctly updated...
  int update_sign(int color) {
    // int new_sign = (hyb_dets->seg_number(color)%2==1 &&
    // (hyb_dets->find(color, 0)).l->dagger )? -1 :1 ;
    int new_sign = (hyb_dets->antiseg_number(color) % 2 == 1) ? -1 : 1;
    int sign_ratio = new_sign * this->sign[color];
    this->sign[color] = new_sign;
    return sign_ratio;
  }

private:
  //-------------------------------------

  /* compute the overlap of a segment
   * \param s: (anti)segment (can start with either type of operator)
   * overlap[i] is overlap of this segment/antisegment with the color i (in
   * particular, <0 for antisegment) \note: the operators are already in the
   * config \note: supposes that the right_occupations are correct!!
   */
  void compute_overlaps(segment const &s);
  double compute_length(int c);

public:
  void check_overlap_matrix_from_scratch();

private:
  double dynamical_contribution(colored_const_iterator cit);
  // compute contribution to trace of a segment from dynamical interactions
  double dynamical_contribution(segment const &s);

  // print
  friend std::ostream &operator<<(std::ostream &out, trace_c_ops const &t);
  friend void h5_write(h5::group fg, std::string subgroup_name,
                       trace_c_ops const &t);
};

} // namespace triqs_ctseg
