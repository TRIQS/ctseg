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
#include "./move_remove_spin_segment.hpp"
namespace triqs_ctseg {
move_remove_spin_segment::move_remove_spin_segment(
    qmc_parameters *params_, configuration *config_,
    triqs::mc_tools::random_generator &RND_)
    : params(params_), config(config_), RND(RND_){};

//--------------------------------------------

double move_remove_spin_segment::attempt() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout
        << "\n =================== ATTEMPT REMOVE SPIN ================ \n"
        << std::endl;
#endif
  tried_insert = false;
  if (config->boson_lines.size() == 0)
    return 0.0;
  int index = RND(config->boson_lines.size());
  auto cseg =
      config->boson_lines[index]; // pick up two spin ops linked by bosonic line

  // bosonic line is not oriented: must determine the 'inside'
  bool no_ops_left_of_plus =
      config->ops_map().cyclic_left(decolor(cseg.plus.sup())) ==
          decolor(cseg.minus.inf()) ||
      config->ops_map().cyclic_left(decolor(cseg.plus.sup())) ==
          decolor(cseg.minus.sup()); // then, remove from plus to minus on the
                                     // left of plus
  bool no_ops_right_of_plus =
      config->ops_map().cyclic_right(decolor(cseg.plus.inf())) ==
          decolor(cseg.minus.inf()) ||
      config->ops_map().cyclic_right(decolor(cseg.plus.inf())) ==
          decolor(cseg.minus.sup()); // then, remove from plus to minus on the
                                     // right of plus

  if (!no_ops_left_of_plus && !no_ops_right_of_plus)
    return 0.0;

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "trying to remove spin segment " << double(cseg.plus.c->tau)
              << ", " << double(cseg.minus.c->tau) << "(index = " << index
              << ")" << std::endl;
  }
#endif
  tried_insert = true;
  double jperp_line_ratio;
  std::tie(jperp_line_ratio, cseg) = config->boson_lines.remove(index);
  auto space_around_seg = (no_ops_right_of_plus)
                              ? config->ops_map().space_around_decolor(
                                    cseg.plus.sup(), cseg.minus.inf())
                              : config->ops_map().space_around_decolor(
                                    cseg.minus.sup(), cseg.plus.inf());

  if (config->ops_map().total_op_number() == 4)
    space_around_seg = double(params->beta); // if last 4 ops...

  double ln_trace_ratio1, ln_trace_ratio2, sign_ratio;
#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "Space_around = " << space_around_seg << std::endl;
    std::cerr << "Removing seg_up" << std::endl;
  }
#endif
  seg_up = no_ops_right_of_plus ? segment{cseg.plus.c, cseg.minus.a}
                                : segment{cseg.minus.a, cseg.plus.c};
  seg_down = no_ops_right_of_plus ? segment{cseg.plus.a, cseg.minus.c}
                                  : segment{cseg.minus.c, cseg.plus.a};

  if (config->ops_map().total_op_number() == 4) {
    double u = RND(1.0);
    seg_up = (u < .5) ? segment{cseg.plus.c, cseg.minus.a}
                      : segment{cseg.minus.a, cseg.plus.c};
    seg_down = (u < .5) ? segment{cseg.plus.a, cseg.minus.c}
                        : segment{cseg.minus.c, cseg.plus.a};
  }
  std::swap(seg_up, seg_down);
  std::tie(ln_trace_ratio1, sign_ratio, seg_desc_up) =
      config->trace.remove_segment(seg_up);

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "Removing seg_down" << std::endl;
  }
#endif
  ln_trace_ratio2 = config->trace.try_remove_segment(seg_down);
  int n_lines_before = config->boson_lines.size() + 1;
  double trace_ratio = std::exp(ln_trace_ratio1 + ln_trace_ratio2);
  double prop_ratio =
      (config->ops_map().total_op_number() == 2)
          ? 2. / (params->beta * params->beta)
          : 2 * n_lines_before / (params->beta * space_around_seg);

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "***RATIOS: TR= " << trace_ratio
              << "|DET= " << jperp_line_ratio << "|PROP = " << prop_ratio
              << " | SIGN " << sign_ratio << std::endl;
  }
#endif

  return trace_ratio * jperp_line_ratio * prop_ratio * sign_ratio;
}

//----------------------------------------------
double move_remove_spin_segment::accept() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "- - - - - ====> ACCEPT - - - - - - - - - - -" << std::endl;
#endif
  config->id++;
  bool seg_up_l_dagger = !seg_down.l->dagger;
  double sign_ratio = config->trace.complete_remove_segment().first;
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
#endif
  // return 1.;
  return sign_ratio;
}

//--------------------------------------------
void move_remove_spin_segment::reject() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "- - - - - ====> REJECT - - - - - - - - - - -" << std::endl;
#endif
  config->id++;
  if (tried_insert) {
    config->trace.reject_remove_segment();
    std::tie(std::ignore, std::ignore, seg_up) =
        config->trace.insert_segment(seg_desc_up);
    config->boson_lines.add({seg_up, seg_down});
  }
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
#endif
}
} // namespace triqs_ctseg
