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
#include "./move_split_spin_segment2.hpp"
namespace triqs_ctseg {
move_split_spin_segment2::move_split_spin_segment2(
    qmc_parameters *params_, configuration *config_,
    triqs::mc_tools::random_generator &RND_)
    : params(params_), config(config_), RND(RND_){};

//--------------------------------------------

double move_split_spin_segment2::attempt() {
  tried_split = false;
  if (config->boson_lines.size() == 0)
    return 0.0;

  // pick up one bosonic line
  int index = RND(config->boson_lines.size());
  auto cseg = config->boson_lines[index];

  B_up = cseg.minus.a; // to be moved
  A_dn = cseg.plus.a;  // to be moved
  B_up_tau = B_up->tau;
  A_dn_tau = A_dn->tau;
  B_dn = cseg.minus.c; // fixed
  A_up = cseg.plus.c;  // fixed

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "trying to split spin segment " << double(cseg.plus.c->tau)
              << ", " << double(cseg.minus.c->tau) << "(index = " << index
              << ")" << std::endl;
  }
#endif
  auto L_up = config->ops_map().cyclic_left(B_up);
  auto R_up = config->ops_map().cyclic_right(B_up);
  auto L_dn = config->ops_map().cyclic_left(A_dn);
  auto R_dn = config->ops_map().cyclic_right(A_dn);

  auto lmax_up = (L_up->tau == R_up->tau) ? params->tau_seg.get_upper_pt()
                                          : L_up->tau - R_up->tau;
  auto lmax_dn = (L_dn->tau == R_dn->tau) ? params->tau_seg.get_upper_pt()
                                          : L_dn->tau - R_dn->tau;

  auto tau_up = params->tau_seg.get_random_pt(RND, lmax_up) + R_up->tau;
  auto tau_dn = params->tau_seg.get_random_pt(RND, lmax_dn) + R_dn->tau;

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "Shifting B_up(" << double(B_up_tau) << ") to "
              << double(tau_up) << ") and A_dn(" << double(A_dn_tau) << " to "
              << double(tau_dn) << std::endl;
  }
#endif

  tried_split = true;

  // remove bosonic line
  double jperp_line_ratio;
  std::tie(jperp_line_ratio, cseg) = config->boson_lines.remove(index);

  // shift G1p and G2p to O1 and O2
  double ln_trace_ratio_up, ln_trace_ratio_dn;

  std::tie(ln_trace_ratio_up, B_up_new) =
      config->trace.try_shift_operator(B_up, tau_up);
  std::tie(ln_trace_ratio_dn, A_dn_new) =
      config->trace.try_shift_operator(A_dn, tau_dn);

  auto det_ratio_up = config->hyb_dets.add(A_up, B_up_new);
  auto det_ratio_dn = config->hyb_dets.try_add(A_dn_new, B_dn);

  double trace_ratio = std::exp(ln_trace_ratio_up + ln_trace_ratio_dn);
  int n_bosonic_lines_before =
      config->boson_lines.size() + 1; // we actually removed a boso line
  int n_hybops_up_after =
      2 * config->hyb_dets.seg_number(0); // actually inserted...
  int n_hybops_dn_after =
      2 * (config->hyb_dets.seg_number(1) + 1); // just tried...

  double prop_ratio = n_bosonic_lines_before * double(lmax_up) *
                      double(lmax_dn) / (n_hybops_up_after * n_hybops_dn_after);

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "***RATIOS: TR= " << trace_ratio
              << "|JPERP= " << jperp_line_ratio << "|PROP = " << prop_ratio
              << "|DET_up= " << det_ratio_up << "|DET_dn= " << det_ratio_dn
              << std::endl;
  }
#endif

  double s_up = (det_ratio_up > 0.0 ? 1.0 : -1.0);
  double s_dn = (det_ratio_dn > 0.0 ? 1.0 : -1.0);
  double prod =
      trace_ratio * jperp_line_ratio * prop_ratio * det_ratio_up * det_ratio_dn;
  return (std::isfinite(prod) ? prod : s_up * s_dn);
}

//----------------------------------------------
double move_split_spin_segment2::accept() {
  config->id++;
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cerr << "ACCEPT SPLIT \n ======" << config->id << "====== \n"
              << std::endl;
#endif
  config->hyb_dets.complete_add();
  double sign_ratio = config->trace.complete_shift_operator();
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
  config->trace.check_overlap_matrix_from_scratch();
#endif
  return sign_ratio;
}

//--------------------------------------------
void move_split_spin_segment2::reject() {
  config->id++;
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cerr << "REJECT SPLIT \n ======" << config->id << "====== \n"
              << std::endl;
#endif
  if (tried_split) {
    config->hyb_dets.reject_add(); // reject addition of G1, G2
    config->hyb_dets.remove(A_up, B_up_new);

    // shift back
    std::tie(std::ignore, B_up) =
        config->trace.try_shift_operator(B_up_new, B_up_tau);
    std::tie(std::ignore, A_dn) =
        config->trace.try_shift_operator(A_dn_new, A_dn_tau);

    // put back bosonic line
    config->boson_lines.add({segment{A_up, B_up}, segment{A_dn, B_dn}});
  }
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
  config->trace.check_overlap_matrix_from_scratch();
#endif
}
} // namespace triqs_ctseg
