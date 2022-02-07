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
#include "./move_group_into_spin_segment2.hpp"
namespace triqs_ctseg {
move_group_into_spin_segment2::move_group_into_spin_segment2(
    qmc_parameters *params_, configuration *config_,
    triqs::mc_tools::random_generator &RND_)
    : params(params_), config(config_), RND(RND_){};

//--------------------------------------------
double move_group_into_spin_segment2::attempt() {
  tried_group = false;
  if (config->hyb_dets.seg_number(0) == 0 ||
      config->hyb_dets.seg_number(1) == 0)
    return 0.0; // no segment to remove
  int op_index_up =
      RND(config->hyb_dets.seg_number(0) * 2); // pick up an operator
  int op_index_dn =
      RND(config->hyb_dets.seg_number(1) * 2); // pick up an operator

  auto seg_up = config->hyb_dets.find(0, op_index_up);
  auto seg_dn = config->hyb_dets.find(1, op_index_dn);

  // two fixed operators (of dagger type)
  if (!seg_up.l->dagger || !seg_dn.l->dagger)
    return 0.0;
  A_up = seg_up.l;
  B_dn = seg_dn.l;

  auto R_dn = config->ops_map().right_neighbor(A_up->tau, 1);
  auto L_dn = config->ops_map().cyclic_left(R_dn);

  auto R_up = config->ops_map().right_neighbor(B_dn->tau, 0);
  auto L_up = config->ops_map().cyclic_left(R_up);

  // two ops to be moved: must be non-dagger!
  A_dn = L_dn->dagger ? R_dn : L_dn;
  B_up = L_up->dagger ? R_up : L_up;

  if (config->boson_lines.contains(A_dn->tau) ||
      config->boson_lines.contains(B_up->tau))
    return 0.0;

  A_dn_tau = A_dn->tau;
  B_up_tau = B_up->tau;

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "trying to group " << double(A_dn->tau) << " to  "
              << double(A_up->tau) << " and " << double(B_up->tau) << " to "
              << double(B_dn->tau) << std::endl;
  }
#endif

  tried_group = true; // all possible trivial rejections have been ruled out
  L_dn = config->ops_map().cyclic_left(A_dn);
  R_dn = config->ops_map().cyclic_right(A_dn);
  L_up = config->ops_map().cyclic_left(B_up);
  R_up = config->ops_map().cyclic_right(B_up);
  auto lmax_dn = (L_dn->tau == R_dn->tau) ? params->tau_seg.get_upper_pt()
                                          : L_dn->tau - R_dn->tau;
  auto lmax_up = (L_up->tau == R_up->tau) ? params->tau_seg.get_upper_pt()
                                          : L_up->tau - R_up->tau;

  // shift O1 and O2 to G1p and G2p
  double ln_trace_ratio_up, ln_trace_ratio_dn;

  auto det_ratio_up = config->hyb_dets.remove(A_up, B_up);
  auto det_ratio_dn = config->hyb_dets.remove(A_dn, B_dn);

  auto eps = params->tau_seg.get_epsilon();
  std::tie(ln_trace_ratio_up, A_dn_new) =
      config->trace.try_shift_operator(A_dn, A_up->tau + eps);
  std::tie(ln_trace_ratio_dn, B_up_new) =
      config->trace.try_shift_operator(B_up, B_dn->tau + eps);

  double jperp_line_ratio = config->boson_lines.add(
      composite_segment(segment{A_up, B_up_new}, segment{A_dn_new, B_dn}));

  double trace_ratio = std::exp(ln_trace_ratio_up + ln_trace_ratio_dn);
  int n_bosonic_lines_after =
      config->boson_lines.size(); // we actually added a boso line
  int n_hybops_up_before =
      2 * (config->hyb_dets.seg_number(0) + 1); // actually removed
  int n_hybops_dn_before =
      2 * (config->hyb_dets.seg_number(1) + 1); // actually removed

  double prop_ratio =
      (n_hybops_up_before * n_hybops_dn_before) /
      (n_bosonic_lines_after * double(lmax_up) * double(lmax_dn));

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "***RATIOS: TR= " << trace_ratio
              << "|JPERP= " << jperp_line_ratio << "|PROP = " << prop_ratio
              << "|DET_dn= " << det_ratio_dn << "|DET_up= " << det_ratio_up
              << std::endl;
  }
#endif

  double s_up = (det_ratio_up > 0.0 ? 1.0 : -1.0);
  double s_dn = (det_ratio_dn > 0.0 ? 1.0 : -1.0);
  double prod =
      trace_ratio * jperp_line_ratio * prop_ratio * det_ratio_up * det_ratio_dn;
  return std::isfinite(prod) ? prod : s_up * s_dn;
}

//----------------------------------------------
double move_group_into_spin_segment2::accept() {
  config->id++;
  double sign_ratio = 1.0;
  if (tried_group) {
    sign_ratio = config->trace.complete_shift_operator();
  }
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
void move_group_into_spin_segment2::reject() {
  config->id++;
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cerr << "REJECT GROUP \n ======" << config->id << "====== \n"
              << std::endl;
#endif
  if (tried_group) {
    config->boson_lines.pop_back();

    // shift back
    std::tie(std::ignore, A_dn) =
        config->trace.try_shift_operator(A_dn_new, A_dn_tau);
    std::tie(std::ignore, B_up) =
        config->trace.try_shift_operator(B_up_new, B_up_tau);

    // config->hyb_dets.reject_remove(); //reject removal of G1, G2
    config->hyb_dets.add(A_dn, B_dn);
    config->hyb_dets.add(A_up, B_up);
  }
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
  config->trace.check_overlap_matrix_from_scratch();
#endif
}
} // namespace triqs_ctseg
