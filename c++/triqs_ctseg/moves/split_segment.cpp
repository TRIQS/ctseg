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
    prop_seg            = sl[prop_seg_idx];
    splitting_full_line = is_full_line(prop_seg);
    if (splitting_full_line) LOG("Splitting full line.");

    // Select splitting points (tau_left,tau_right)
    tau_t dt1 = tau_t::random(rng, prop_seg.length());
    tau_t dt2 = tau_t::random(rng, prop_seg.length());
    if (dt1 == dt2) {
      LOG("Generated equal times");
      return 0;
    }
    if (dt1 > dt2 and !splitting_full_line)
      std::swap(dt1, dt2); // If splitting a full line, the order of tau_left and tau_right is not fixed
    tau_left             = prop_seg.tau_c - dt1; // dt1 < dt2
    tau_right            = prop_seg.tau_c - dt2;
    auto removed_segment = segment_t{tau_left, tau_right}; // "antisegment" : careful with order of c, cdag
    // Check whether the splitting of a cyclic segment produces a new segment at the beginning of seglist
    // (useful for config update and dets)
    bool segment_overboard = is_cyclic(prop_seg) and !is_cyclic(segment_t{tau_right, prop_seg.tau_cdag});
    // Index of the rightmost of the two produced segments (the one to the left is always at prop_seg_idx)
    right_seg_idx = (splitting_full_line or segment_overboard) ? 0 : prop_seg_idx + 1;

    LOG("Splitting at position {} : adding c at {}, cdag at {}", prop_seg_idx, tau_right, tau_left);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = -wdata.mu(color) * removed_segment.length();
    for (auto c : range(config.n_color())) {
      if (c != color) { ln_trace_ratio -= -wdata.U(color, c) * overlap(config.seglists[c], removed_segment); }
      if (wdata.has_Dt) { ln_trace_ratio += K_overlap(config.seglists[c], tau_right, tau_left, wdata.K, color, c); }
    }
    if (wdata.has_Dt)
      ln_trace_ratio += -real(wdata.K(double(tau_left - tau_right))(color, color)); // Correct double counting
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    auto &D             = wdata.dets[color];
    auto det_c_time     = [&](long i) { return D.get_y(i).first; };
    auto det_cdag_time  = [&](long i) { return D.get_x(i).first; };
    long det_index_c    = lower_bound(det_c_time, D.size(), tau_right);
    long det_index_cdag = lower_bound(det_cdag_time, D.size(), tau_left);

    // We insert tau_cdag as a line (first index) and tau_c as a column (second index).
    auto det_ratio = D.try_insert(det_index_cdag, det_index_c, {tau_left, 0}, {tau_right, 0});

    // ------------  Proposition ratio ------------

    double current_number_segments = sl.size();
    double future_number_intervals = splitting_full_line ? 1 : double(sl.size()) + 1.0;
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

    double initial_sign = config_sign(config, wdata.dets);
    LOG("Initial sign is {}. Initial configuration: {}", initial_sign, config);

    // Update the dets
    wdata.dets[color].complete_operation();

    // Split the segment
    auto &sl = config.seglists[color];
    if (splitting_full_line) {
      auto new_segment = segment_t{tau_right, tau_left};
      sl[prop_seg_idx] = new_segment;
    } else {
      auto new_seg_left  = segment_t{prop_seg.tau_c, tau_left, prop_seg.J_c, false};
      auto new_seg_right = segment_t{tau_right, prop_seg.tau_cdag, false, prop_seg.J_cdag};
      // Update the proposed segment
      sl[prop_seg_idx] = new_seg_left;
      // Insert a new segment
      sl.insert(sl.begin() + right_seg_idx, new_seg_right);
    }

    double final_sign = config_sign(config, wdata.dets);
    double sign_ratio = final_sign / initial_sign;
    LOG("Final sign is {}", final_sign);

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata.dets);
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}. Config: {}.",
                   det_sign, sign_ratio, config);

    LOG("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void split_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
}; // namespace moves
