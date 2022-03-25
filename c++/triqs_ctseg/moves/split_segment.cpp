#include "split_segment.hpp"
#include "../logs.hpp"

namespace moves {

  double split_segment::attempt() {

    LOG("\n =================== ATTEMPT SPLIT ================ \n");

    // ------------ Choice of segment --------------
    // Select color
    color    = rng(wdata.n_color);
    auto &sl = config.seglists[color];
    LOG("Splitting at color {}", color);

    if (sl.empty()) {
      LOG("Nothing to split");
      return 0;
    }

    // Select segment to split
    prop_seg_idx        = rng(sl.size());
    prop_seg            = sl[prop_seg_idx];
    splitting_full_line = is_full_line(prop_seg, fac);
    if (splitting_full_line) LOG("Splitting full line.");

    // Select splitting points (tau_left,tau_right)
    qmc_time_t prop_seg_length = splitting_full_line ? wdata.qmc_beta : prop_seg.tau_c - prop_seg.tau_cdag;
    qmc_time_t dt1             = wdata.make_random_time(rng, prop_seg_length);
    qmc_time_t dt2             = wdata.make_random_time(rng, prop_seg_length);
    if (dt1 == dt2) {
      LOG("Generated equal times");
      return 0;
    }
    if (dt1 > dt2 && !splitting_full_line)
      std::swap(dt1, dt2); // If splitting a full line, the order of tau_left and tau_right is not fixed
    tau_left             = prop_seg.tau_c - dt1; // dt1 < dt2
    tau_right            = prop_seg.tau_c - dt2;
    auto removed_segment = segment_t{tau_left, tau_right}; // "antisegment" : careful with order of c, cdag
    // Check whether the splitting of a cyclic segment produces a new segment at the beginning of seglist
    // (useful for config update and dets)
    segment_overboard = is_cyclic(prop_seg) and !is_cyclic(segment_t{tau_right, prop_seg.tau_cdag});
    // Index of the rightmost of the two produced segments (the one to the left is always at prop_seg_idx)
    right_seg_idx = (splitting_full_line or segment_overboard) ? 0 : prop_seg_idx + 1;

    LOG("Splitting at position {} : adding c at {}, cdag at {}", prop_seg_idx, tau_right, tau_left);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = -wdata.mu(color) * removed_segment.length();
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio -= -wdata.U(color, c) * overlap(config.seglists[c], removed_segment, fac);
        if (wdata.has_Dt)
          ln_trace_ratio -= K_overlap(config.seglists[c], removed_segment, slice_target_to_scalar(wdata.K, color, c));
      }
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    /* We insert tau_cdag as a line (first index) and tau_c as a column (second index). The index always corresponds to the
     segment the tau_c/tau_cdag belongs to. Here, the cdag is always inserted at the position of the segment we are splitting.
     The c insertion position depends on whether we are splitting a full line and whether the right segment is "overboard", computed 
     in right_seg_idx. */
    auto det_ratio = wdata.dets[color].try_insert(prop_seg_idx, right_seg_idx, {tau_left, 0}, {tau_right, 0});
    // additional sign due to "roll" is not computed at this stage, cf accept

    // ------------  Proposition ratio ------------

    double current_number_segments = sl.size();
    double future_number_intervals = splitting_full_line ? 1 : double(sl.size()) + 1.0;
    double prop_ratio = (current_number_segments * prop_seg_length * prop_seg_length / (splitting_full_line ? 1 : 2))
       / (future_number_intervals);

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double split_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Update the dets
    wdata.dets[color].complete_operation();

    auto &sl          = config.seglists[color];
    double sign_ratio = 1;
    // Split the segment and compute the sign ratio
    if (splitting_full_line) {
      auto new_segment = segment_t{tau_right, tau_left};
      if (is_cyclic(new_segment)) sign_ratio = -1;
      sl[prop_seg_idx] = new_segment;
    } else {
      auto new_seg_left  = segment_t{prop_seg.tau_c, tau_left};
      auto new_seg_right = segment_t{tau_right, prop_seg.tau_cdag};
      if (is_cyclic(prop_seg)) {
        bool kept_number_cyclic = is_cyclic(new_seg_left) or is_cyclic(new_seg_right);
        // If we are splitting a cyclic segment and both new segments are not cyclic, get a - sign
        if (not kept_number_cyclic) sign_ratio *= -1;
        if (segment_overboard) sign_ratio *= wdata.dets[color].roll_matrix(det_t::Down);
      }

      LOG("Sign ratio is {}", sign_ratio);
      // Update the proposed segment
      sl[prop_seg_idx] = new_seg_left;
      // Insert a new segment
      sl.insert(sl.begin() + right_seg_idx, new_seg_right);
    }

    // Check invariant
#ifdef CHECK_INVARIANTS
    check_invariant(config, wdata.dets);
#endif
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}. Config: {}.",
                   det_sign, sign_ratio, config);

    SPDLOG_TRACE("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void split_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
}; // namespace moves
