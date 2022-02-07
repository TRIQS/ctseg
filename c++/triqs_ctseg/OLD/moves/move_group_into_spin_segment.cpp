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
#include "./move_group_into_spin_segment.hpp"
namespace triqs_ctseg {
move_group_into_spin_segment::move_group_into_spin_segment(
    qmc_parameters *params_, configuration *config_,
    triqs::mc_tools::random_generator &RND_)
    : params(params_), config(config_), RND(RND_){};

//--------------------------------------------
/// find point to group with
/**
 * traverse from L to R (to the right) looking for operator G such that G is
 * - closest to O
 * - on opposite line
 * - with opposite dagger
 * - neighbor of O
 * @return {found_G, G}
 */
std::pair<bool, colored_const_iterator>
move_group_into_spin_segment::find_point_to_group_with(
    colored_const_iterator O, colored_const_iterator L,
    colored_const_iterator R) {

  auto R_dec = config->ops_map().cyclic_right(decolor(O));
  if (R_dec->tau != R->tau && R_dec->dagger != O->dagger)
    return {true, config->ops_map().find_colored(R_dec->color, R_dec->tau)};
  auto L_dec = config->ops_map().cyclic_left(decolor(O));
  if (L_dec->tau != L->tau && L_dec->dagger != O->dagger)
    return {true, config->ops_map().find_colored(L_dec->color, L_dec->tau)};
  return {false, R};
}

//--------------------------------------------
double move_group_into_spin_segment::attempt() {
  tried_group = false;
  int color = RND(params->n_color); // pick up color
  if (config->hyb_dets.seg_number(color) == 0)
    return 0.0; // no segment to remove
  int op_index_1 =
      RND(config->hyb_dets.seg_number(color) * 2); // pick up an operator
  int op_index_2 =
      RND(config->hyb_dets.seg_number(color) * 2); // pick up an operator
#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "[grouping] color = " << color
              << ", op_index_1 = " << op_index_1 << std::endl;
  }
#endif

  if (op_index_1 == op_index_2)
    return 0.0;

  auto seg1 = config->hyb_dets.find(color, op_index_1);
  auto seg2 = config->hyb_dets.find(color, op_index_2);
  O1 = seg1.l;
  O2 = seg2.l;
  if (O1->dagger == O2->dagger)
    return 0.0; // O1 and O2 must be of opposite type

  auto L1 = config->ops_map().cyclic_left(O1);
  auto R1 = config->ops_map().cyclic_right(O1);
  auto L2 = config->ops_map().cyclic_left(O2);
  auto R2 = config->ops_map().cyclic_right(O2);

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "[group] selected times" << double(O1->tau) << " and  "
              << double(O2->tau) << " with n.n (" << double(L1->tau) << ", "
              << double(R1->tau) << ") and (" << double(L2->tau) << ", "
              << double(R2->tau) << ")" << std::endl;
  }
#endif

  bool found_G1;
  colored_const_iterator G1;
  std::tie(found_G1, G1) = find_point_to_group_with(O1, L1, R1);
  if (!found_G1)
    return 0.0;
  bool found_G2;
  colored_const_iterator G2;
  std::tie(found_G2, G2) = find_point_to_group_with(O2, L2, R2);
  if (!found_G2)
    return 0.0;

  //  if(config->hyb_dets.seg_number(1-color)==0) return 0.0;
  //  auto L1_dec = config->ops_map().cyclic_left(decolor(O1));
  //  while(L1_dec->color!=(1-color)){ L1_dec =
  //  config->ops_map().cyclic_left(L1_dec);} auto L1 = config->ops_map().

  bool O1_O2_neighbors = (O1 == config->ops_map().cyclic_left(O2)) ||
                         (O1 == config->ops_map().cyclic_right(O2));
  bool G1_G2_neighbors = (G1 == config->ops_map().cyclic_left(G2)) ||
                         (G1 == config->ops_map().cyclic_right(G2));

  if (O1_O2_neighbors || G1_G2_neighbors)
    return 0.0; // exclude case when G1/G2 **OR** O1/O2 are adjacent
  // if both O1/O2 and G1/G2 have intercalated ops, the split move must choose
  // among both lines otherwise, no choice of line is to be made
  bool no_choice_of_line = (O1_O2_neighbors != G1_G2_neighbors);

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "trying to group " << double(O1->tau) << " to  "
              << double(G1->tau) << " and " << double(O2->tau) << " to  "
              << double(G2->tau) << std::endl;
  }
#endif
  tried_group = true; // all possible trivial rejections have been ruled out
  auto lmax1 = L1->tau - R1->tau;
  auto lmax2 = L2->tau - R2->tau;

  // shift O1 and O2 to G1p and G2p
  double ln_trace_ratio1, ln_trace_ratio2;

  O1_tau = O1->tau;
  O2_tau = O2->tau;

  auto det_ratio_1 =
      config->hyb_dets.remove(O1, O2); // actually removes one hyb_segment //
                                       // could do try here: not same line!
  auto det_ratio_2 = config->hyb_dets.try_remove(G1, G2);

  // careful to epsilon...
  auto eps = params->tau_seg.get_epsilon();
  std::tie(ln_trace_ratio1, G1p) =
      config->trace.try_shift_operator(O1, G1->tau + eps);
  std::tie(ln_trace_ratio2, G2p) =
      config->trace.try_shift_operator(O2, G2->tau + eps);

  double jperp_line_ratio = config->boson_lines.add(
      composite_segment(segment{G1, G2}, segment{G1p, G2p}));

  double trace_ratio = std::exp(ln_trace_ratio1 + ln_trace_ratio2);
  int n_bosonic_lines_after =
      config->boson_lines.size(); // we actually added a boso line
  int n_hybops_before = 2 * (config->hyb_dets.seg_number(color) +
                             1); // we haven't touched the dets
  double prop_ratio =
      (no_choice_of_line)
          ? (config->n_colors() / 2 * n_hybops_before * n_hybops_before) /
                (n_bosonic_lines_after * double(lmax1) * double(lmax2))
          : (config->n_colors() / 2 * n_hybops_before * n_hybops_before) /
                (2 * n_bosonic_lines_after * double(lmax1) * double(lmax2));

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "***RATIOS: TR= " << trace_ratio
              << "|JPERP= " << jperp_line_ratio << "|PROP = " << prop_ratio
              << "|DET1= " << det_ratio_1 << "|DET2= " << det_ratio_2
              << std::endl;
  }
#endif

  double s_1 = (det_ratio_1 > 0.0 ? 1.0 : -1.0);
  double s_2 = (det_ratio_2 > 0.0 ? 1.0 : -1.0);
  double prod =
      trace_ratio * jperp_line_ratio * prop_ratio * det_ratio_1 * det_ratio_2;
  return std::isfinite(prod) ? prod : s_1 * s_2;
}

//----------------------------------------------
double move_group_into_spin_segment::accept() {
  config->id++;
  config->hyb_dets.complete_remove();
  double sign_ratio = config->trace.complete_shift_operator();
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cerr << "ACCEPT GROUP (sign = " << sign_ratio
              << ") \n ======" << config->id << "====== \n"
              << std::endl;
  if (config->print_condition())
    config->print_to_h5();
  config->trace.check_overlap_matrix_from_scratch();
#endif
  return sign_ratio;
}

//--------------------------------------------
void move_group_into_spin_segment::reject() {
  config->id++;
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cerr << "REJECT GROUP \n ======" << config->id << "====== \n"
              << std::endl;
#endif
  if (tried_group) {
    config->boson_lines.pop_back();

    // shift back
    std::tie(std::ignore, O1) = config->trace.try_shift_operator(G1p, O1_tau);
    std::tie(std::ignore, O2) = config->trace.try_shift_operator(G2p, O2_tau);

    config->hyb_dets.reject_remove(); // reject removal of G1, G2
    config->hyb_dets.add(O1, O2);     // add back
  }
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
  config->trace.check_overlap_matrix_from_scratch();
#endif
}
} // namespace triqs_ctseg
