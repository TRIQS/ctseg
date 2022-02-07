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
#include "./move_split_spin_segment.hpp"
namespace triqs_ctseg {
move_split_spin_segment::move_split_spin_segment(
    qmc_parameters *params_, configuration *config_,
    triqs::mc_tools::random_generator &RND_)
    : params(params_), config(config_), RND(RND_){};

//--------------------------------------------

double move_split_spin_segment::attempt() {
  tried_insert = false;
  if (config->boson_lines.size() == 0)
    return 0.0;
  if (config->ops_map().total_op_number() == 4)
    return 0.0; // easily handled by removal update

  // pick up one bosonic line
  int index = RND(config->boson_lines.size());
  auto cseg = config->boson_lines[index];

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "trying to split spin segment " << double(cseg.plus.c->tau)
              << ", " << double(cseg.minus.c->tau) << "(index = " << index
              << ")" << std::endl;
  }
#endif

  // choose line [up/down] where ops will be moved around
  bool up_has_no_ops_left =
      config->ops_map().cyclic_left(cseg.plus.c) == cseg.minus.a;
  bool up_has_no_ops_right =
      config->ops_map().cyclic_right(cseg.plus.c) == cseg.minus.a;
  bool down_has_no_ops_right =
      config->ops_map().cyclic_left(cseg.minus.c) == cseg.plus.a;
  bool down_has_no_ops_left =
      config->ops_map().cyclic_right(cseg.minus.c) == cseg.plus.a;

  if (up_has_no_ops_left || down_has_no_ops_left || up_has_no_ops_right ||
      down_has_no_ops_right)
    return 0.0;
  bool no_choice_of_line = false;
  int color = bool(RND(2));

  seg_not_shifted = (color == 0) ? segment{cseg.plus.a, cseg.minus.c}
                                 : segment{cseg.plus.c, cseg.minus.a};

  // pick up ops G1p and G2p and compute space around them
  G1p = (color == 0) ? cseg.plus.c : cseg.plus.a;
  G2p = (color == 0) ? cseg.minus.a : cseg.minus.c;
  auto G1 = (color == 0) ? cseg.plus.a : cseg.plus.c;
  auto G2 = (color == 0) ? cseg.minus.c : cseg.minus.a;
#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "G1p = " << double(G1p->tau) << std::endl;
    std::cerr << "G2p = " << double(G2p->tau) << std::endl;
  }
#endif
  G1p_tau = G1p->tau;
  G2p_tau = G2p->tau;

  auto R1 = config->ops_map().cyclic_right(G1p);
  auto L1 = config->ops_map().cyclic_left(G1p);
  auto R2 = config->ops_map().cyclic_right(G2p);
  auto L2 = config->ops_map().cyclic_left(G2p);

  auto lmax1 = L1->tau - R1->tau;
  auto lmax2 = L2->tau - R2->tau;

  auto tau_1 = params->tau_seg.get_random_pt(RND, lmax1) + R1->tau;
  auto tau_2 = params->tau_seg.get_random_pt(RND, lmax2) + R2->tau;
#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "Shifting G1p(" << double(G1p_tau) << ") and G2p("
              << double(G2p_tau) << ") to " << double(tau_1) << " and "
              << double(tau_2) << std::endl;
  }
#endif

  tried_insert = true;

  // remove bosonic line
  double jperp_line_ratio;
  std::tie(jperp_line_ratio, cseg) = config->boson_lines.remove(index);

  // shift G1p and G2p to O1 and O2
  double ln_trace_ratio1;
  double ln_trace_ratio2;

  std::tie(ln_trace_ratio1, O1) = config->trace.try_shift_operator(G1p, tau_1);
  std::tie(ln_trace_ratio2, O2) = config->trace.try_shift_operator(G2p, tau_2);

  auto det_ratio_1 = config->hyb_dets.add(O1, O2);
  auto det_ratio_2 = config->hyb_dets.try_add(G1, G2);

  double trace_ratio = std::exp(ln_trace_ratio1 + ln_trace_ratio2);
  int n_bosonic_lines_before =
      config->boson_lines.size() + 1; // we actually removed a boso line
  int n_hybops_after =
      2 * config->hyb_dets.seg_number(color) + 4; // we have inserted 4 hyb_ops

  double prop_ratio =
      (no_choice_of_line)
          ? n_bosonic_lines_before * double(lmax1) * double(lmax2) /
                (config->n_colors() / 2 * n_hybops_after * n_hybops_after)
          : 2 * n_bosonic_lines_before * double(lmax1) * double(lmax2) /
                (config->n_colors() / 2 * n_hybops_after * n_hybops_after);

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
  return (std::isfinite(prod) ? prod : s_1 * s_2);
}

//----------------------------------------------
double move_split_spin_segment::accept() {
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
void move_split_spin_segment::reject() {
  config->id++;
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cerr << "REJECT SPLIT \n ======" << config->id << "====== \n"
              << std::endl;
#endif
  if (tried_insert) {
    config->hyb_dets.reject_add(); // reject addition of G1, G2
    config->hyb_dets.remove(O1, O2);

    // shift back
    std::tie(std::ignore, G1p) = config->trace.try_shift_operator(O1, G1p_tau);
    std::tie(std::ignore, G2p) = config->trace.try_shift_operator(O2, G2p_tau);

    // put back bosonic line
    config->boson_lines.add({seg_not_shifted, segment{G1p, G2p}});
  }
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
  config->trace.check_overlap_matrix_from_scratch();
#endif
}
} // namespace triqs_ctseg
