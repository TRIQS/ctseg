#include "split_segment.hpp"
#include "../logs.hpp"

namespace moves {

  double split_segment::attempt() {

    LOG("\n =================== ATTEMPT SPLIT ================ \n");

    // ------------ Choice of segment --------------
    // Select color
    color    = rng(config.n_color());
    auto &sl = config.seglists[color];
    LOG("Splitting at color {}", color);

    if (sl.empty()) {
      LOG("Nothing to split");
      return 0;
    }

    // Select segment to split
    prop_seg_idx        = rng(sl.size());
    auto prop_seg       = sl[prop_seg_idx];
    splitting_full_line = is_full_line(prop_seg);
    if (splitting_full_line) LOG("Splitting full line.");

    // Select splitting points (tau_left,tau_right)
    auto dt1 = tau_t::random(rng, prop_seg.length());
    auto dt2 = tau_t::random(rng, prop_seg.length());
    if (dt1 == dt2) {
      LOG("Generated equal times");
      return 0;
    }
    // If splitting a full line, the order of tau_left and tau_right is not fixed
    if (dt1 > dt2 and not splitting_full_line) std::swap(dt1, dt2);

    tau_left  = prop_seg.tau_c - dt1; // dt1 < dt2
    tau_right = prop_seg.tau_c - dt2;

    LOG("Splitting at position {} : adding c at {}, cdag at {}", prop_seg_idx, tau_right, tau_left);

    // ------------  Trace ratio  -------------

    auto removed_segment = segment_t{tau_left, tau_right}; // "antisegment" : careful with order of c, cdag

    double ln_trace_ratio = -wdata.mu(color) * removed_segment.length();
    for (auto c : range(config.n_color())) {
      if (c != color) { ln_trace_ratio -= -wdata.U(color, c) * overlap(config.seglists[c], removed_segment); }
      if (wdata.has_Dt) { ln_trace_ratio += K_overlap(config.seglists[c], tau_right, tau_left, wdata.K, color, c); }
    }
    if (wdata.has_Dt)
      ln_trace_ratio += -real(wdata.K(double(tau_left - tau_right))(color, color)); // Correct double counting
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    auto &bl     = wdata.block_number[color];
    auto &bl_idx = wdata.index_in_block[color];
    auto &D      = wdata.dets[bl];
    if (wdata.offdiag_delta) {
      if (cdag_in_det(tau_left, D) or c_in_det(tau_right, D)) {
        LOG("One of the proposed times already exists in another line of the same block. Rejecting.");
        return 0;
      }
    }
    auto det_ratio = D.try_insert(det_lower_bound_x(D, tau_left),  //
                                  det_lower_bound_y(D, tau_right), //
                                  {tau_left, bl_idx}, {tau_right, bl_idx});

    // ------------  Proposition ratio ------------

    double current_number_segments = sl.size();
    double future_number_intervals = splitting_full_line ? 1 : sl.size() + 1;
    // T direct 1/  # segment  1/ len(prop_seg) ^2 * (2 iif !full line)
    // T inverse : 1/ # interval
    double prop_ratio =
       (current_number_segments * prop_seg.length() * prop_seg.length() / (splitting_full_line ? 1 : 2))
       / (future_number_intervals);

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    double prod = trace_ratio * det_ratio * prop_ratio;

    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double split_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    double initial_sign = config_sign(wdata);
    LOG("Initial sign is {}. Initial configuration: {}", initial_sign, config);

    // Update the dets
    wdata.dets[wdata.block_number[color]].complete_operation();

    // Split the segment
    auto &sl = config.seglists[color];
    if (splitting_full_line) {
      sl[prop_seg_idx] = segment_t{tau_right, tau_left};
    } else {
      auto prop_seg      = sl[prop_seg_idx];
      auto new_seg_left  = segment_t{prop_seg.tau_c, tau_left, prop_seg.J_c, false};
      auto new_seg_right = segment_t{tau_right, prop_seg.tau_cdag, false, prop_seg.J_cdag};
      // Update the proposed segment
      sl[prop_seg_idx] = new_seg_left;
      // Insert the new seg at the right
      // if the splitting of a cyclic segment produces a non cyclic one, i.e. at the front, we need to insert it at the
      // front of the list
      bool insert_at_front = is_cyclic(prop_seg) and not is_cyclic(new_seg_right);
      sl.insert(sl.begin() + (insert_at_front ? 0 : prop_seg_idx + 1), new_seg_right);
    }

    double final_sign = config_sign(wdata);
    double sign_ratio = final_sign / initial_sign;
    LOG("Final sign is {}", final_sign);

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata);

    if (sign_ratio * det_sign == -1.0) wdata.minus_sign = true;

    LOG("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void split_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[wdata.block_number[color]].reject_last_try();
  }
}; // namespace moves
