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
#include "./move_insert_spin_segment.hpp"
namespace triqs_ctseg {

move_insert_spin_segment::move_insert_spin_segment(
    const qmc_parameters *params_, configuration *config_,
    triqs::mc_tools::random_generator &RND_)
    : params(params_), config(config_), RND(RND_), special_case(false){};

//----------------------------------------------
/// from a time, color, compute (length to the left neighbour, length to the
/// right neighbour, iterator to the right neighbour
// if no operator, return {beta, beta, end(color)}
std::tuple<qmc_time_t, qmc_time_t, colored_const_iterator>
move_insert_spin_segment::max_lens(qmc_time_t const &tau, int color) {
  if (config->ops_map().op_number(color) == 0)
    return std::make_tuple(params->tau_seg.get_upper_pt(),
                           params->tau_seg.get_upper_pt(),
                           config->ops_map().end(color));
  auto R = config->ops_map().right_neighbor(tau, color);
  if (R == config->ops_map().end(color))
    return std::make_tuple(params->tau_seg.get_upper_pt(),
                           params->tau_seg.get_upper_pt(), R); // no op
  auto L = config->ops_map().cyclic_left(R);
  return std::make_tuple(L->tau - tau, tau - R->tau, R);
}

//----------------------------------------------
double move_insert_spin_segment::attempt() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "\n =================== ATTEMPT INSERT SPIN================ \n"
              << std::endl;
#endif

  // ----------  Selection of two times -------------------
  auto tau_1 = params->tau_seg.get_random_pt(
      RND); // choose tau between 0 and beta: position of operator

  // find maximal lengths to te right and left on up and do
  // line+right_neighbors.
  colored_const_iterator R_up, R_do;
  qmc_time_t l_max_up_l, l_max_up_r, l_max_do_l, l_max_do_r;
  std::tie(l_max_up_l, l_max_up_r, R_up) = max_lens(
      tau_1,
      0); // compute l_max and iterator to right neighbor of tau_1, line up
  std::tie(l_max_do_l, l_max_do_r, R_do) = max_lens(
      tau_1,
      1); // compute l_max and iterator to right neighbor of tau_1, line do

  // determine interval where tau_2 can be inserted
  auto l_max_l = std::min(l_max_up_l, l_max_do_l);
  auto l_max_r = std::min(l_max_up_r, l_max_do_r);
  qmc_time_t lmax = (config->ops_map().seg_number(0) == 0 &&
                     config->ops_map().seg_number(1) == 0)
                        ? params->tau_seg.get_upper_pt()
                        : l_max_l + l_max_r;
  qmc_time_t delta_length = params->tau_seg.get_random_pt(RND, lmax);
  auto eps = params->tau_seg.get_epsilon();
  qmc_time_t tau_2 = tau_1 + l_max_l - delta_length; // careful to ordering...

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "Trying to insert spin at tau_1 = " << tau_1
              << " and tau_2 = " << tau_2 << " (lmax = " << double(lmax) << ")"
              << std::endl;
  }
#endif

  // rejection if trying to insert tau_1 where both lines have same occupation
  bool up_occupied = (config->ops_map().seg_number(0) > 0)
                         ? R_up->dagger
                         : config->trace.full_lines(0);
  bool do_occupied = (config->ops_map().seg_number(1) > 0)
                         ? R_do->dagger
                         : config->trace.full_lines(1);
  if (up_occupied == do_occupied) {
    special_case = true;
    return 0.0;
  }

  // swap: none if no ops at all, otherwise fixed by right neighbor R
  auto R = (config->ops_map().seg_number(0) == 0)
               ? R_do
               : R_up; // if no ops on line 0, line 1 will fix need_swap
  bool need_swap = (config->ops_map().total_op_number() == 0)
                       ? false
                       : (on_same_side(tau_1, tau_2, R->tau) ? !(tau_1 > tau_2)
                                                             : (tau_1 > tau_2));
  if (need_swap)
    std::swap(tau_1, tau_2);

  auto sd1 = segment_desc{tau_1, tau_2, 0};
  auto sd2 = segment_desc{tau_1 - eps, tau_2 + eps, 1};
#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "sd1 = " << sd1.tl << ", " << sd1.tr << std::endl;
    std::cerr << "sd2 = " << sd2.tl << ", " << sd2.tr << std::endl;
  }
#endif
  // insertion
  double ln_trace_ratio1, ln_trace_ratio2, sign_ratio1;
  try {
    std::tie(ln_trace_ratio1, sign_ratio1, seg1) =
        config->trace.insert_segment(sd1);
  } catch (insertion_error const &e) {
    special_case = true;
    std::cerr << "insert error -- move insert spin" << std::endl;
    return 0.0;
  }

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    bool is_plus = seg1.l->dagger;
    if (is_plus)
      std::cerr << "Spin segment" << std::endl;
    else
      std::cerr << "Spin anti-segment" << std::endl;
    std::cerr << "Inserted segment " << double(sd1.tl) << "," << double(sd1.tr)
              << " on line 0 (seg1.l->dagger=" << seg1.l->dagger << ")"
              << std::endl;
  }
#endif
  segment seg2;
  try {
    std::tie(ln_trace_ratio2, seg2) = config->trace.try_insert_segment(sd2);
  } catch (insertion_error const &e) {
    special_case = true;
    config->trace.remove_segment(seg1);
    return 0.0;
  };

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "Inserted segment " << double(sd2.tl) << "," << double(sd2.tr)
              << std::endl;
  }
#endif
  int n_lines_after = config->boson_lines.size() + 1;
  double trace_ratio = std::exp(ln_trace_ratio1 + ln_trace_ratio2);
#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "ln1=" << ln_trace_ratio1 << std::endl;
    std::cerr << "ln2=" << ln_trace_ratio2 << std::endl;
  }
#endif
  double prop_ratio = (config->ops_map().total_op_number() == 2)
                          ? params->beta * params->beta / 2
                          : params->beta * lmax / (2 * n_lines_after);
  double boson_ratio = config->boson_lines.add(composite_segment(seg1, seg2));

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "RATIOS: TR = " << trace_ratio << "| PROP = " << prop_ratio
              << "| BLINE = " << boson_ratio
              << "| SIGN_RATIO_1 = " << sign_ratio1 << std::endl;
  }
#endif

  return trace_ratio * prop_ratio * boson_ratio * sign_ratio1;
}

//----------------------------------------------
double move_insert_spin_segment::accept() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "- - - - - ====> ACCEPT - - - - - - - - - - -" << std::endl;
#endif
  config->id++;
  double sign_ratio = config->trace.complete_insert_segment();
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
  config->trace.check_overlap_matrix_from_scratch();
#endif
  special_case = false;
  // return 1.;
  return sign_ratio;
}

//----------------------------------------------
void move_insert_spin_segment::reject() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "- - - - - ====> REJECT - - - - - - - - - - -" << std::endl;
#endif
  config->id++;
  if (!special_case) {
    config->boson_lines.pop_back();
    config->trace.reject_insert_segment();
    config->trace.remove_segment(seg1);
  }
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
  config->trace.check_overlap_matrix_from_scratch();
#endif
  special_case = false;
}
} // namespace triqs_ctseg
