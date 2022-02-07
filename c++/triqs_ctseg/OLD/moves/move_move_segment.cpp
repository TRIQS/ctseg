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
#include "./move_move_segment.hpp"
namespace triqs_ctseg {

/*
 * move which takes one (anti)segment from one line to another line
 * Move is decomposed into an removal and insert move
 */

double move_move_segment::attempt() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "\n =================== ATTEMPT MOVE ================ \n"
              << std::endl;
#endif
  move_performed = false;

  // --------  Selection of the times and color for the segment to be moved
  // ------------

  // choose two lines: one original color, one final color
  c1 = RND(params->n_color);
  c2 = RND(params->n_color);
  if (c2 == c1)
    return 0.0; // c2 = (c1+1)%(params->n_color);
  // check if there are operators on the selected lines
  if (config->hyb_dets.seg_number(c1) == 0)
    return 0.0;
  auto seg1 =
      config->hyb_dets.find(c1, RND(config->hyb_dets.seg_number(c1) *
                                    2)); // find corresponding (anti)segment
  // if overlap on c2 > 0 : cannot move ->return 0.0
  if (config->ops_map().op_count_in_between(seg1.l, seg1.r, c2) > 0)
    return 0.0; // if the c2 has operators where we want to move seg1, reject.

  auto R2 = config->ops_map().right_neighbor(seg1.r->tau, c2);
  if (config->ops_map().seg_number(c2) > 0 && R2->dagger != seg1.l->dagger)
    return 0.0; // trying to move (anti)segment where there is already a
                // (anti)segment
  // if (config->ops_map().seg_number(c2)==0 && (config->trace.full_lines(c2) !=
  // seg1.l->dagger)) return 0.0; //prevent moving (anti)segment to (empty) full
  // line
  if (config->hyb_dets.seg_number(c2) == 0 &&
      (config->trace.full_lines(c2) != seg1.l->dagger))
    return 0.0; // prevent moving (anti)segment to (empty) full line

#ifdef CTSEG_DEBUG
  std::cout << "Trying to move segment tau_l = " << seg1.l->tau
            << ", tau_r = " << seg1.r->tau << " of line " << c1 << " to line "
            << c2 << std::endl;
  std::cout << " config before the move is " << *config << std::endl;
#endif

  const int n1 = config->hyb_dets.seg_number(c1),
            n2 = config->hyb_dets.seg_number(c2);
  // --------  Removing the first segment, try to insert the second and get
  // trace and det ratio ---------------
  double det_ratio, ln_trace_ratio1, ln_trace_ratio2;

  diff_block = (config->color_to_block_and_inner_index(c1).first !=
                config->color_to_block_and_inner_index(c2).first);

  // removing the first segment
  if (diff_block)
    det_ratio = config->hyb_dets.try_remove2(
        seg1.l, seg1.r); // actually removes seg1 from hyb_opmap
  else
    det_ratio = config->hyb_dets.remove(seg1.l, seg1.r);

  std::tie(ln_trace_ratio1, sign_ratio, seg_d1) =
      config->trace.remove_segment(seg1);

  // try inserting the second segment
  segment seg2;
  std::tie(ln_trace_ratio2, seg2) =
      config->trace.try_insert_segment(segment_desc{seg_d1.tl, seg_d1.tr, c2});
  det_ratio *= config->hyb_dets.try_add(seg2.l, seg2.r);

  double trace_ratio = std::exp(ln_trace_ratio1 + ln_trace_ratio2);

  // -------- Proposition ratio and finish ---------------

  double prop_ratio = double(n1) / (n2 + 1.);
  double s = (det_ratio > 0.0 ? 1.0 : -1.0);
  double prod = trace_ratio * det_ratio * prop_ratio;
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cerr << "***RATIOS: TR= " << trace_ratio << "|DET= " << det_ratio
              << "|PROP = " << prop_ratio << std::endl;
#endif
  move_performed = true;
  return (std::isfinite(prod) ? prod : s);
}

//--------------------------------------------------

double move_move_segment::accept() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "- - - - - ====> ACCEPT - - - - - - - - - - -" << std::endl;
#endif
  config->id++;
  if (diff_block)
    config->hyb_dets.complete_remove2();
  // else config->hyb_dets.complete_remove();//does not do anything
  config->hyb_dets.complete_add(); // adds seg2 to hyb_opmap
  sign_ratio *= config->trace.complete_insert_segment();

#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
#endif
  return sign_ratio;
}

//--------------------------------------------------

void move_move_segment::reject() {
  config->id++;
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "- - - - - ====> REJECT - - - - - - - - - - -" << std::endl;
  std::cerr << "config: " << *config << std::endl;
#endif
  if (!move_performed)
    return;
  config->hyb_dets.reject_add();
  config->trace.reject_insert_segment();

  // put back the first segment
  segment seg;
  std::tie(std::ignore, std::ignore, seg) = config->trace.insert_segment(
      seg_d1); // inserts back, but sign is wrong because segment has actually
               // been removed from hyb_opmap
  if (diff_block)
    config->hyb_dets.reject_remove2(seg.l, seg.r);
  else
    config->hyb_dets.add(seg.l, seg.r);
  config->trace.update_sign(c1);
  // config->hyb_dets.add(seg.l,seg.r);
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
  std::cerr << "config: " << *config << std::endl;
#endif
}

} // namespace triqs_ctseg
