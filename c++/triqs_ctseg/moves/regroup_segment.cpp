#include "regroup_segment.hpp"
#include "../logs.hpp"

namespace moves {

  double regroup_segment::attempt() {

    LOG("\n =================== ATTEMPT REGROUP ================ \n");

    // ------------ Choice of segment --------------
    // Select color
    color    = rng(config.n_color());
    auto &sl = config.seglists[color];
    LOG("Regrouping at color {}", color);

    // If no segments nothing to regroup
    if (sl.empty()) {
      LOG("Color is empty");
      return 0;
    }

    // Select pair of segments (or cyclic segment) to regroup
    making_full_line = sl.size() == 1;
    if (making_full_line) {
      if (is_full_line(sl[0])) {
        LOG("Segment is full line");
        return 0; // If segment is full line nothing to regroup
      }
      left_seg_idx  = 0;
      right_seg_idx = 0;
    } else {
      left_seg_idx  = rng(sl.size());
      right_seg_idx = (left_seg_idx == sl.size() - 1) ? 0 : left_seg_idx + 1;
    }
    left_seg  = sl[left_seg_idx];
    right_seg = sl[right_seg_idx];

    if (left_seg.J_cdag or right_seg.J_c) {
      LOG("At least one of the operators has spin line attached, cannot regroup.");
      return 0;
    }

    LOG("Regroup at positions {} and {}: removing c at {}, cdag at {}", left_seg_idx, right_seg_idx, right_seg.tau_c,
        left_seg.tau_cdag);

    // ------------  Trace ratio  -------------

    auto inserted_seg = segment_t{left_seg.tau_cdag, right_seg.tau_c}; // "antisegment" : careful with order of c, cdag
    double ln_trace_ratio = wdata.mu(color) * inserted_seg.length();

    for (auto c : range(config.n_color())) {
      if (c != color) { ln_trace_ratio += -wdata.U(color, c) * overlap(config.seglists[c], inserted_seg); }
      if (wdata.has_Dt) {
        ln_trace_ratio -= K_overlap(config.seglists[c], right_seg.tau_c, left_seg.tau_cdag, wdata.K, color, c);
      }
    }
    if (wdata.has_Dt)
      ln_trace_ratio -=
         real(wdata.K(double(right_seg.tau_c - left_seg.tau_cdag))(color, color)); // Correct double counting

    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    // We remove a cdag (first index) from the left segment and a c (second index) from the right segment.
    auto bl        = wdata.block_number[color];
    auto &D        = wdata.dets[bl];
    auto det_ratio = D.try_remove(det_lower_bound_x(D, left_seg.tau_cdag), //
                                  det_lower_bound_y(D, right_seg.tau_c));

    // ------------  Proposition ratio ------------

    double future_number_segments   = making_full_line ? 1 : sl.size() - 1;
    double current_number_intervals = sl.size();
    // Length of future segment
    auto new_seg_len = (making_full_line ? tau_t::beta() : left_seg.tau_c - right_seg.tau_cdag);

    // T direct  = 1/current_number_intervals
    // T inverse = 1/ future_number_segments / len_of_new_seg ^2 * (2 iif !full line)
    double prop_ratio =
       current_number_intervals / (future_number_segments * new_seg_len * new_seg_len / (making_full_line ? 1 : 2));

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    double prod = trace_ratio * det_ratio * prop_ratio;

    return (std::isfinite(prod)) ? prod : det_sign;
  }

  //--------------------------------------------------

  double regroup_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    double initial_sign = trace_sign(wdata);
    LOG("Initial sign is {}. Initial configuration: {}", initial_sign, config);

    // Update the dets
    wdata.dets[wdata.block_number[color]].complete_operation();

    // Regroup segments
    auto &sl = config.seglists[color];
    if (making_full_line) {
      sl[left_seg_idx] = segment_t{tau_t::beta(), tau_t::zero()};
    } else {
      // Update the left segment
      sl[left_seg_idx] = segment_t{left_seg.tau_c, right_seg.tau_cdag, left_seg.J_c, right_seg.J_cdag};
      // Remove the right segment
      sl.erase(sl.begin() + right_seg_idx);
    }

    double final_sign = trace_sign(wdata);
    double sign_ratio = final_sign / initial_sign;
    LOG("Final sign is {}", final_sign);

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata);

    if (sign_ratio * det_sign == -1.0) wdata.minus_sign = true;
    LOG("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void regroup_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[wdata.block_number[color]].reject_last_try();
  }
} // namespace moves
